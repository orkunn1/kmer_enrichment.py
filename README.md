# kmer_enrichment.py
# kmer-enrichment

Compare k-mer frequencies between two FASTA sequence sets (target vs. background) and test for statistical enrichment using Fisher's exact test in Python.

## Why This Tool?

When comparing a target sequence set (e.g., edited, perturbed, or candidate regions) against a background/control set, raw counts alone do not show whether a short sequence (`k-mer`) is overrepresented by chance or statistically enriched. This single-script tool calculates k-mer counts, frequencies, enrichment ratios, and one-sided Fisher's exact `p_value`s in one command.

## Requirements

- `pandas`
- `scipy`

Install with:

```bash
pip install pandas scipy
