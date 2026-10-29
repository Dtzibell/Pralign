from typing import Any
import json
import requests as rq
from ctypes import ArgumentError
import polars as pl
import parasail as ps
import pathlib
import argparse
import re
import os
import tomllib
from .api import get_tm_batch_onetaq_stdbuf, get_tm_batch_phusion_hfbuf

# DNA scoring matrix
MATRIX = ps.matrix_create("ACGT", 1, -1)
COMP = str.maketrans("ACGT", "TGCA")
VALID_QUERY = re.compile("^[A|C|G|T]+$")
LABGURU = "https://mftuebingen.labguru.com/api/v1/primers"


def reverse_and_complement(query: str) -> str:
    """Makes a complement and reverses the given query sequence.
    :param query: A query DNA sequence.
    """
    return query.translate(COMP)[::-1]


# TODO: would it not make sense to make gap_open infinity? Would avoid filtering down the line
def align(query: str, seq: str, gap_open: int = 1) -> dict[str, int | str]:
    """Performs Smith-Waterman on the given query and sequence.
    :param query: DNA sequence 1.
    :param seq: DNA sequence 2.
    :param gap_open: penalty for opening and extending a gap.
    """
    result = ps.sw_trace(query, seq, gap_open, gap_open, MATRIX)
    aligned_primer = result.traceback.ref
    return {"alignment score": result.score, "aligned portion": aligned_primer}


def load_dotenv(path=".env"):
    """Minimal .env loader: KEY=value lines, '#' comments, optional quotes.
    Real environment variables take precedence over the file."""
    if not os.path.isfile(path):
        return
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key = key.removeprefix("export ").strip()
            value = value.strip().strip('"').strip("'")
            os.environ.setdefault(key, value)


def update_db(silent_errors: bool) -> None:
    print("Updating the database")
    load_dotenv()
    token = os.environ.get("API_KEY")
    if token is None:
        raise ArgumentError(
            "Could not find API_KEY in .env, hence updating the database is not possible"
        )
    response = True
    counter = 1
    df = pl.DataFrame(
        {"primer name": [], "vial number": [], "sequence": []},
        schema={
            "primer name": pl.String,
            "vial number": pl.String,
            "sequence": pl.String,
        },
    )
    while response:
        print(f"Retrieving page {counter}")
        response = rq.get(
            f"{LABGURU}",
            params={
                "token": token,
                "content_type": "application/json",
                "page": counter,
            },
        ).json()
        counter += 1
        titles = []
        sequences = []
        vial_nos = []
        for primer in response:
            try:
                sequences.append(primer["sequences"][0]["seq"])
                titles.append(primer["title"])
                vial_nos.append(primer["Vial Number"])
            except IndexError as e:
                if not silent_errors:
                    if hasattr(e, "message"):
                        print(
                            f"{e.message}. Likely the primer {primer['title']} contains no sequence"
                        )
                    else:
                        print(
                            f"{e}. Likely the primer {primer['title']} contains no sequence"
                        )
        titles = pl.Series("primer name", titles)
        vial_nos = pl.Series("vial number", vial_nos)
        sequences = pl.Series("sequence", sequences)
        df = pl.concat([df, pl.DataFrame([titles, vial_nos, sequences])])
    df.write_parquet("database.parquet", compression="zstd", compression_level=22)
    print(f"Updated the database. Total primers: {df.height}")


def main() -> None:
    with open("config.toml", "rb") as f:
        toml = tomllib.load(f)
    argp = argparse.ArgumentParser(
        prog="Pralign",
        description="Find a primer sequence that matches an arbitrarily long query DNA sequence",
    )
    argp.add_argument("query")
    argp.add_argument(
        "--update",
        "-u",
        help="Whether to update the existing database",
        action="store_true",
        default=toml["database"]["update"],
    )
    argp.add_argument(
        "--reverse",
        "-r",
        help="Whether to reverse and make a complement of the query strand",
        action="store_true",
        default=toml["query"]["reverse"],
    )
    argp.add_argument(
        "--top",
        "-n",
        help="The top N primers to output into a csv file",
        default=toml["output"]["primers"]["top"],
    )
    argp.add_argument(
        "--longest",
        "-l",
        help="The length of the longest primer of interest",
        default=toml["filter"]["primers"]["longest"],
    )
    argp.add_argument(
        "--shortest",
        "-s",
        help="The length of the shortest primer of interest",
        default=toml["filter"]["primers"]["shortest"],
    )
    argp.add_argument(
        "--output",
        "-o",
        help='The name of the output file, without the extension. E.g. "alignments"',
        default=toml["output"]["filename"],
    )
    namespace = argp.parse_args()
    run(namespace, toml)


def parse_query(query: str, reverse: bool) -> str:
    query = query.strip().upper()
    if VALID_QUERY.fullmatch(query) is None:
        raise ArgumentError(
            f"query {query} is not valid. Only characters A, C, T, and G are allowed"
        )
    if reverse:
        query = reverse_and_complement(query)
    return query


def parse_shortest(shortest: str) -> int:
    shortest: int = int(shortest)
    if shortest < 2:
        raise ArgumentError(
            f'The argument "shortest" must be at least 2. Given {shortest}.'
        )
    return shortest


def run(namespace: argparse.Namespace, toml: dict[str, Any]) -> None:
    if namespace.update:
        update_db(toml["runtime"]["silent_errors"])

    query = parse_query(namespace.query, namespace.reverse)
    top: int = int(namespace.top)
    longest: int = int(namespace.longest)
    shortest = parse_shortest(namespace.shortest)
    output: str = namespace.output.strip()

    print("Aligning primers")
    db: pathlib.Path = pathlib.Path("database.parquet")
    df: pl.DataFrame = (
        pl.scan_parquet(db)
        .with_columns(pl.col("sequence").str.len_bytes().alias("sequence length"))
        .filter(
            (pl.col("sequence length") < longest)
            & (pl.col("sequence length") > shortest)
        )  # some sequences are question marks? Not relevant, but must set lower bound to at least 2
        .with_columns(
            pl.col("sequence")
            .map_elements(
                lambda seq: align(query, seq.upper()),
                return_dtype=pl.Struct(
                    {"alignment score": pl.Int64, "aligned portion": pl.String}
                ),
            )
            .alias("alignment")
        )  # performs the alignment
        .unnest("alignment")  # unnests the struct into two separate columns
        .filter(
            pl.col("aligned portion")
            .str.contains("-")
            .not_()
            .and_(pl.col("alignment score") > 15)
        )  # filters out gapped and insufficient alignments. TODO: introduce a toml for params? Overkill?
        .sort(by="alignment score", descending=True)
        .with_columns(
            (pl.col("alignment score") / pl.col("sequence length")).alias(
                "aligned proportion"
            )
        )  # Computes the aligned proportion of the primer. Alignment score = alignment length because match score is 1.
        .sort(by="aligned proportion", descending=True)
        .head(top)
        .collect()
    )
    onetaq_tm = pl.Series(
        "OneTaq Tm",
        [
            dct["tm1"]
            for dct in get_tm_batch_onetaq_stdbuf(df["aligned portion"].to_list())
        ],
    )
    phusion_tm = pl.Series(
        "Phusion Tm",
        [
            dct["tm1"]
            for dct in get_tm_batch_phusion_hfbuf(df["aligned portion"].to_list())
        ],
    )
    df = df.with_columns(onetaq_tm, phusion_tm).drop("aligned portion")
    df.write_csv(f"{output}.csv", float_precision=3)
    print(f"Alignment output written to {output}.csv")
