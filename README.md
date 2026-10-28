# Pralign
LabGuru primer database searcher. \
Performs local alignment on a given query sequence against a primer database to find primers that could possibly fit the sequence. \
The script is linked to the [NEB Tm calculator](https://tmcalculator.neb.com/#!/main) via their [Tm calculation API](https://tmapi.neb.com/), so the aligned primers immediately come with their respective Tm's.

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
  * `--longest` or `-l`: The length of the longest primer of interest.
  * `--shortest` or `-s`: The length of the shortest primer of interest.
  * `--output` or `-o`: The name of the output file, without the extension. E.g. `alignments`

And a flag:

  * `--reverse` or `-r`: Whether to reverse input sequence and make the complement of the input strand the query. (Since copying the complementary strand is annoying in SnapGene) 

Example argument usage: 
```bash
uv run pralign ACTGACTGACTGACTGACTGACTGACTGACTG \
--top 20 --longest 18 --shortest 40 \
--output alignments --reverse
``` 
Coincidentally, the provided values are the defaults. \
If you ever forget these arguments and flags, you can always review them with `uv run pralign -h`. 

## Notes on usage
The primers are sorted by the aligned proportion of their sequence, so the most likely fits are at the top. Gapped alignments are excluded from the results.

Finding primers for multiple sequences is relatively easy with a bash script:
```bash
#!/usr/bin/bash
uv run pralign ACTG -o meaningfulAlignmentName1
uv run pralign ACTG -n 10 -r \
--output meaningfulReversedAlignmentName2WithN10
```
or a powershell script:
```powershell
uv run pralign ACTG -o meaningfulAlignmentName1
uv run pralign ACTG -n 10 -r \
-o meaningfulReversedAlignmentName2WithN10
```
The script is not well-tested. Reach out to me or your favourite AI assistant if something's unclear.

Unfortunately SnapGene primer import is only available with the paid version and the documentation is poor, so I did not make the files importable.

You can save multiple databases in the same directory. The `.xlsx` file chosen as the database is the one that is lexicographically sorted to be the last. So `2026-10-28_Primer.xlsx` will be chosen over `2026-09-28_Primer.xlsx` if that is of any relevance
