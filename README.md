# Pralign
LabGuru primer database searcher. Performs local alignment on a given query sequence against a primer database to find primers that could possibly fit the sequence.

## Prerequisites
Officially Python 3.14 if you know what you are doing, but I suggest installing [uv](https://docs.astral.sh/uv/#installation). It's typically quick and painless.

## How to use

  * Clone the repository with `git clone https://github.com/Dtzibell/Pralign.git`
  * Export the primer database from LabGuru: Inventory $\rightarrow$ Primers $\rightarrow$ ... (Meatball menu) $\rightarrow$ Export Sequences.
  * Move the `.xslx` file to the project's top level directory (where the `pyproject.toml` is located).
  * Run `uv sync` once to download necessary python libraries
  * Run `uv run pralign ACTG`, where `ACTG` may be any DNA query that you want to find a primer for.

The programme takes further optional arguments:

  * `--reverse` or `-r`: reverse and make the complement of the input as the query. (Since its annoying to make 5' $\rightarrow$ 3' copies in SnapGene)

