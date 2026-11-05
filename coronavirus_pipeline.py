#!/usr/bin/env python3
"""Unified translation of the user's six uploaded Perl/R scripts.

Requires: NCBI BLAST+ (tblastn), samtools, and Biopython (pip install biopython).

The main pipeline reproduces the intended order of pipe.pl.  The extra
sequence-extraction and BLAST/FASTA-filter routines embedded in blastadjust(1).pl
are available as separate subcommands.
"""
import argparse
import csv
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path

FIELDS = ('qseqid','sseqid','pident','length','mismatch','gapopen',
          'qstart','qend','sstart','send','evalue','bitscore')


def run(cmd, stdout=None):
    print('+ ' + ' '.join(map(str, cmd)), file=sys.stderr)
    return subprocess.run([str(x) for x in cmd], stdout=stdout, check=True)


def blast(query, db, dest, max_targets):
    with dest.open('w') as out:
        run(['tblastn', '-query', query, '-db', db, '-outfmt', '6',
             '-max_target_seqs', str(max_targets)], stdout=out)


def read_blast(path):
    with Path(path).open(newline='') as handle:
        for n, row in enumerate(csv.reader(handle, delimiter='\t'), 1):
            if not row or row[0].startswith('#'):
                continue
            if len(row) != 12:
                raise ValueError(f'{path}:{n}: expected 12 BLAST outfmt 6 columns, got {len(row)}')
            yield dict(zip(FIELDS, row))


def write_hits(path, hits):
    # Clean, headerless outfmt-6 table. Unlike R write.table, no row names.
    with Path(path).open('w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        for hit in hits:
            writer.writerow([hit[k] for k in FIELDS])


def select_hits(blast_file, matches, mismatches):
    """blastselect.r: highest bitscore per subject; qstart==1 split.

    Ties retain the first BLAST hit for that subject.
    """
    best = OrderedDict()
    for hit in read_blast(blast_file):
        subject = hit['sseqid']
        if subject not in best or float(hit['bitscore']) > float(best[subject]['bitscore']):
            best[subject] = hit
    good = [h for h in best.values() if int(h['qstart']) == 1]
    bad = [h for h in best.values() if int(h['qstart']) != 1]
    write_hits(matches, good)
    write_hits(mismatches, bad)
    return len(good), len(bad)


def ensure_index(fasta):
    if not Path(str(fasta) + '.fai').exists():
        run(['samtools', 'faidx', fasta])


def fasta_reference_ids(fasta):
    with Path(str(fasta) + '.fai').open() as f:
        return {line.split('\t', 1)[0]: int(line.split('\t')[1]) for line in f if line.strip()}


def extract_fasta(fasta, hits_file, output, protein_length=None, reverse_complement=False):
    """Extract BLAST subject intervals, optionally extending to target AA length.

    Original blastadjust.pl used mismatches as alignment length, and query ID
    rather than subject ID. Both are corrected here. Coordinates are 1-based.
    Reverse hits are extended toward lower coordinates; optionally reverse-
    complemented to coding orientation.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    ensure_index(fasta)
    lengths = fasta_reference_ids(fasta)
    records = []
    for hit in read_blast(hits_file):
        subject = hit['sseqid']
        if subject not in lengths:
            raise ValueError(f'Subject {subject!r} absent from FASTA index for {fasta}')
        sstart, send = int(hit['sstart']), int(hit['send'])
        is_reverse = sstart > send
        lo, hi = min(sstart, send), max(sstart, send)
        if protein_length is not None:
            # Match the original calculation: (expected AA - aligned AA) * 3.
            extension = (protein_length - int(hit['length'])) * 3
            if extension < 0:
                print(f'Warning: {subject}: alignment longer than expected protein; no extension', file=sys.stderr)
                extension = 0
            if is_reverse:
                lo -= extension
            else:
                hi += extension
        lo = max(1, lo)
        hi = min(lengths[subject], hi)
        if lo > hi:
            print(f'Warning: skipping empty region {subject}:{lo}-{hi}', file=sys.stderr)
            continue
        region = f'{subject}:{lo}-{hi}'
        result = subprocess.run(['samtools', 'faidx', str(fasta), region],
                                text=True, capture_output=True, check=True)
        from io import StringIO
        rec = SeqIO.read(StringIO(result.stdout), 'fasta')
        if reverse_complement and is_reverse:
            rec.seq = rec.seq.reverse_complement()
            rec.description = rec.id
        records.append(rec)
    with Path(output).open('w') as out:
        SeqIO.write(records, out, 'fasta')
    return len(records)


def simplify_fasta(src, dst):
    """Literal Fastasimplifier.pl and Fastasimplifier2.pl header operation.

    Every header is truncated at its FIRST period; sequence lines unchanged.
    """
    with Path(src).open() as inp, Path(dst).open('w') as out:
        for line in inp:
            if line.startswith('>'):
                out.write(line.split('.', 1)[0].rstrip('\r\n') + '\n')
            else:
                out.write(line)


def genbank_annotations(genbank, fasta, output):
    """parsegb.pl: output FASTA header, GenBank DEFINITION, /organism.

    Match on the full header (minus '>') as the original Perl intended.
    Also accept the accession's versionless form because simplification
    removes .version. Unmatched headers produce blank fields.
    """
    from Bio import SeqIO
    lookup = {}
    for record in SeqIO.parse(str(genbank), 'genbank'):
        definition = record.description
        organism = record.annotations.get('organism', '')
        ids = {record.id, record.id.split('.')[0]}
        for accession in record.annotations.get('accessions', []):
            ids.add(accession)
            ids.add(accession.split('.')[0])
        for key in ids:
            lookup[key] = (definition, organism)
    with Path(fasta).open() as inp, Path(output).open('w') as out:
        for line in inp:
            if line.startswith('>'):
                name = line[1:].strip()
                definition, organism = lookup.get(name, ('', ''))
                out.write(f'{name}\t{definition}\t{organism}\n')


def filter_fasta(fasta, blast_file, output, min_identity=25, min_length=1000):
    """Optional SequenceGP code from blastadjust(1).pl.

    Select subject IDs where identity > threshold AND length > threshold.
    The original also dumped every input FASTA record before the matches;
    this version writes only selected sequences.
    """
    from Bio import SeqIO
    records = {r.id: r for r in SeqIO.parse(str(fasta), 'fasta')}
    selected = OrderedDict()
    for hit in read_blast(blast_file):
        if float(hit['pident']) > min_identity and int(hit['length']) > min_length:
            subject = hit['sseqid']
            if subject in records:
                selected[subject] = records[subject]
            else:
                print(f'Warning: {subject} not found in {fasta}', file=sys.stderr)
    with Path(output).open('w') as out:
        SeqIO.write(selected.values(), out, 'fasta')
    return len(selected)


def pipeline(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    blast_file = outdir / 'blasttable.txt'
    matches = outdir / 'blastmatches.txt'
    mismatches = outdir / 'blastmismatches.txt'
    raw = outdir / 'out.fa'
    final1 = outdir / 'Finalout.fasta'
    final2 = outdir / 'Finalout2.fasta'
    annotations = outdir / 'annotations.txt'
    if args.skip_blast:
        if not blast_file.is_file():
            raise FileNotFoundError(f'--skip-blast requires {blast_file}')
    else:
        blast(args.query, args.database, blast_file, args.max_targets)
    good, bad = select_hits(blast_file, matches, mismatches)
    print(f'Selected {good} qstart=1 hits; {bad} other best hits', file=sys.stderr)
    count = extract_fasta(args.reference, matches, raw, args.length,
                          reverse_complement=args.reverse_complement)
    simplify_fasta(raw, final1)
    genbank_annotations(args.genbank, final1, annotations)
    simplify_fasta(final1, final2)
    print(f'Extracted {count} sequences. Outputs: {outdir}', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('run', help='Full pipeline corresponding to pipe.pl')
    p.add_argument('query', help='Protein query FASTA')
    p.add_argument('database', help='tblastn nucleotide BLAST database prefix')
    p.add_argument('genbank', help='GenBank file containing matching accessions')
    p.add_argument('length', type=int, help='Expected protein length in amino acids')
    p.add_argument('--reference', default='Allcoronas2.fasta', help='Nucleotide FASTA matching BLAST subject IDs')
    p.add_argument('--outdir', default='.', help='Output directory')
    p.add_argument('--max-targets', type=int, default=500)
    p.add_argument('--skip-blast', action='store_true', help='Reuse outdir/blasttable.txt')
    p.add_argument('--reverse-complement', action='store_true', help='Orient reverse-strand hits to coding strand')
    p.set_defaults(func=pipeline)
    p = sub.add_parser('select', help='Run Python equivalent of blastselect.r')
    p.add_argument('blast_table')
    p.add_argument('--matches', default='blastmatches.txt')
    p.add_argument('--mismatches', default='blastmismatches.txt')
    p.set_defaults(func=lambda a: print(select_hits(a.blast_table, a.matches, a.mismatches)))
    p = sub.add_parser('extract', help='Exact or extended samtools extraction')
    p.add_argument('blast_table')
    p.add_argument('reference', help='FASTA file; e.g. Newseqs2.fasta')
    p.add_argument('output')
    p.add_argument('--length', type=int, help='Expected protein length (AA); omit for exact extraction')
    p.add_argument('--reverse-complement', action='store_true')
    p.set_defaults(func=lambda a: print(f'Extracted {extract_fasta(a.reference, a.blast_table, a.output, a.length, a.reverse_complement)} records'))
    p = sub.add_parser('filter', help='Optional SequenceGP identity/length FASTA filter')
    p.add_argument('fasta')
    p.add_argument('blast_table')
    p.add_argument('output')
    p.add_argument('--min-identity', type=float, default=25)
    p.add_argument('--min-length', type=int, default=1000)
    p.set_defaults(func=lambda a: print(f'Selected {filter_fasta(a.fasta, a.blast_table, a.output, a.min_identity, a.min_length)} records'))
    p = sub.add_parser('simplify', help='Original FASTA header simplifier')
    p.add_argument('input')
    p.add_argument('output')
    p.set_defaults(func=lambda a: simplify_fasta(a.input, a.output))
    p = sub.add_parser('annotate', help='GenBank accession/definition/organism mapping')
    p.add_argument('genbank')
    p.add_argument('fasta')
    p.add_argument('output')
    p.set_defaults(func=lambda a: genbank_annotations(a.genbank, a.fasta, a.output))
    args = parser.parse_args()
    try:
        args.func(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'Error: {exc}\n')


if __name__ == '__main__':
    main()