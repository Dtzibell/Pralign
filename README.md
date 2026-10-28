# Pralign
LabGuru primer database searcher. Performs local alignment on a given query sequence against a primer database to find primers that could possibly fit the sequence.

## Prerequisites
Officially Python 3.14 if you know what you are doing, but I suggest installing [uv](https://docs.astral.sh/uv/#installation). It's typically quick and painless.

## How to use

  * Clone the repository with `git clone https://github.com/Dtzibell/Pralign.git`
  * Export the primer database from LabGuru: Inventory $\rightarrow$ Primers $\rightarrow$ ... (Meatball menu) $\rightarrow$ Export Sequences. Rename the primer database to `primers_<DATE_OF_EXPORT>.xlsx`
  * Move the `.xslx` file to the project's top level directory (where the `pyproject.toml` is located).
  * Run `uv sync` once to download necessary python libraries
  * Run `uv run pralign ACTG`, where `ACTG` may be any DNA query that you want to find a primer for.

The programme takes further optional arguments:

  * `--reverse` or `-r`: reverse and make the complement of the input as the query. (Since its annoying to make 5' $\rightarrow$ 3' copies in SnapGene)

The output is a file named `alignments.csv` which contains the name, the vial number, the sequence, the length of the sequence (`seq_len`) and the score (`align_score`) of the top 20 best scoring primers.
## Notes on usage
I recommend selecting primers whos `align_score` and `seq_len` are at the worst very close to each other.

Most of the variables are currently hard-coded, notably the lengths of the primers that are are filtered. If it is of interest to find primers that are larger than around 20 bp, I recommend to play around with the numbers in `seq/pralign/__init__.py`. Reach out to me or your favourite AI assistant if something is unclear.

You can save multiple databases in the same directory. The `.xlsx` file chosen as the database is the one that is lexicographically sorted to be the last. So `primers_2026_10_28.xlsx` will be chosen over `primers_2026_10_27.xslx` if that is of any relevance
