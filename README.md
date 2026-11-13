# Pralign
Have you ever wondered if your primer database already has a primer that suits a segment that you are trying to amplify? Pralign answers exactly this question. \
Pralign performs local alignment on a given query sequence against a primer database to find primers that could possibly fit the sequence. 

Pralign is integrated with the LabGuru API, so updating the primer database is quick.\
Additionally, the script is linked to the [NEB Tm calculator](https://tmcalculator.neb.com/#!/main) via their [Tm calculation API](https://tmapi.neb.com/), so the aligned primers immediately come with their respective Tm's.

## Prerequisites
I suggest installing [uv](https://docs.astral.sh/uv/#installation). It's typically quick and painless.

## Quick start

  * Clone the repository with `git clone https://github.com/Dtzibell/Pralign.git`
  * Create an `.env` file. Add a line `API_KEY=<insert your API key here>` to the file. To generate an API key, follow LabGuru's [tutorial](https://labguru-support.cenevo.com/hc/en-us/articles/22786313894940-How-to-generate-a-personal-token).
  * Run `uv run pralign`. 
  
This opens a GUI, with which you can perform database search. Update your database, enter your sequence into the bottom left and press "Align"!

## Command line interface
The programme also support a command line interface. A minimal command is `uv run pralign <query-sequence>`.

The programme takes additional optional arguments:

  * `--top` or `-n`: The top N primers to output into a csv file.
  * `--longest` or `-l`: The length of the longest primer of interest.
  * `--shortest` or `-s`: The length of the shortest primer of interest.
  * `--output` or `-o`: The name of the output file, without the extension. E.g. `alignments`

And flags:

  * `--reverse` or `-r`: Whether to reverse input sequence and make the complement of the input strand the query. (Since copying the complementary strand is annoying in SnapGene) 
  * `--update` or `-u`: Whether to update the database during this run.

Example argument usage: 

```bash
uv run pralign ACTGACTGACTGACTGACTGACTGACTGACTG \
--top 20 --longest 17 --shortest 40 \
--output alignments --update --reverse
``` 

Coincidentally, the provided values are the defaults. \
If you ever forget these arguments and flags, you can always review them with `uv run pralign -h`. 

If you would like to change the default parameters, you can modify `config.toml`. The format is quite intuitive, but refer to [TOML documentation](https://toml.io/en/) if you crash due to invalid values.

## Notes on usage
The primers are sorted by the aligned proportion of their sequence, so the most likely fits are at the top. The programme allows at most one gap, which is then replaced by an N in the $T_{m}$ calculation.

Finding primers for multiple sequences is relatively easy with a script:
```bash
uv run pralign ACTG -o meaningfulAlignmentName1
uv run pralign ACTG -n 10 -r \
--output meaningfulReversedAlignmentName2WithN10

The script is not well-tested. Reach out to me or your favourite AI assistant if something's unclear.

Unfortunately SnapGene primer import is only available with the paid version and the documentation is poor, so I did not make the files importable.
