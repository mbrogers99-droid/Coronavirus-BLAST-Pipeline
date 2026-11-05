# Coronavirus-BLAST-Pipeline

A Python workflow for extracting coronavirus nucleotide regions corresponding to protein-query `tblastn` hits, selecting best hits, simplifying FASTA identifiers, and matching GenBank accession metadata.

This consolidates the author's original Perl and R scripts into one Python command-line program. It has been used by the author to reproduce results from the previous workflow; independent testing across other datasets has not yet been performed.

## Requirements

- Python 3.9 or newer
- [Biopython](https://biopython.org/) (`pip install -r requirements.txt`)
- NCBI BLAST+ (`tblastn`)
- samtools (`samtools faidx`)

Both external executables must be available on `PATH`. Your BLAST nucleotide database must already be formatted with `makeblastdb`, and the reference FASTA must contain the corresponding BLAST subject IDs.

## Usage

```bash
python coronavirus_pipeline.py run \
    RBD_query.fasta \
    CoronavirusDB \
    Coronavirus.gb \
    222 \
    --reference Allcoronas2.fasta \
    --outdir results
```

The four positional arguments are the protein query FASTA, BLAST nucleotide database prefix, GenBank input, and expected protein length in amino acids. **Replace `222` with the correct length for your query.**

### Main outputs

| Output | Description |
| --- | --- |
| `blasttable.txt` | BLAST outfmt 6 hits |
| `blastmatches.txt` | Highest-bit-score hit per subject where query start is 1 |
| `blastmismatches.txt` | Highest-bit-score hit per subject where query start is not 1 |
| `out.fa` | Extracted nucleotide sequences |
| `Finalout.fasta` | FASTA with simplified headers |
| `Finalout2.fasta` | Second simplified FASTA output |
| `annotations.txt` | Tab-separated accession, definition, and organism annotations |

### Other commands

```bash
python coronavirus_pipeline.py --help
python coronavirus_pipeline.py select --help
python coronavirus_pipeline.py extract --help
python coronavirus_pipeline.py filter --help
python coronavirus_pipeline.py simplify --help
python coronavirus_pipeline.py annotate --help
```

`--reverse-complement` is available for orienting reverse-strand hits. Without it, extracted sequence orientation follows the reference FASTA.

## Notes and limitations

- Best hits are chosen by BLAST bit score per subject sequence, retaining the first hit on ties.
- The main pipeline processes hits with `qstart == 1` and writes other best hits separately.
- The extraction extension depends on the expected protein length; verify that this matches your input protein query.
- `annotations.txt` contains GenBank accession metadata, **not** a per-gene feature table.
- Large FASTA, BLAST database, and result files are excluded from version control by `.gitignore`.

## License

MIT License. See [LICENSE](LICENSE).