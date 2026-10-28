# Pralign
LabGuru primer database searcher. Performs local alignment on a given query sequence against a primer database to find primers that could possibly fit the sequence.

## Prerequisites
Officially Python 3.14 if you know what you are doing, but I suggest installing [uv](https://docs.astral.sh/uv/#installation). It's typically quick and painless.

## Quick start

  * Clone the repository with `git clone https://github.com/Dtzibell/Pralign.git`
  * Export the primer database from LabGuru: Inventory $\rightarrow$ Primers $\rightarrow$ ... (Meatball menu) $\rightarrow$ Export Sequences.
  * Move the `.xslx` file to the project's top level directory (where the `pyproject.toml` is located).
  * Run `uv sync` once to download necessary python libraries
  * Run `uv run pralign ACTG`, where `ACTG` may be any DNA query that you want to find a primer for.

The output is a file named `alignments.csv` which contains the name, the vial number, the sequence, the length of the sequence (`seq_len`) and the score (`align_score`) of the top 20 best scoring primers.

## Optional arguments 

The programme takes additional optional arguments:

  * `--top` or `-n`: The top N primers to output into a csv file.
  * `--longest` or `-l`: The of the longest primer of interest.
  * `--shortest` or `-s`: The of the shortest primer of interest.

And a flag:

  * `--reverse` or `-r`: Whether to reverse and make a complement of the query strand reverse and make the complement of the input as the query (Since copying the complementary strand is annoying in SnapGene). 

Example argument usage: `uv run pralign ACTGACTGACTGACTGACTGACTGACTGACTG --top 20 --longest 18 --shortest 40 --reverse`. Coincidentally, the provided values are the defaults. \
If you ever forget these arguments and flags, you can always review them with `uv run pralign -h`. 

## Notes on usage
I recommend selecting primers whos `align_score` and `seq_len` are at the worst very close to each other.

The script is not well-tested. Reach out to me or your favourite AI assistant if something's unclear.

You can save multiple databases in the same directory. The `.xlsx` file chosen as the database is the one that is lexicographically sorted to be the last. So `2026-10-28_Primer.xlsx` will be chosen over `2026-09-28_Primer.xlsx` if that is of any relevance
