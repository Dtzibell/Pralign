from ctypes import ArgumentError
import polars as pl
import parasail as ps
import pathlib
import argparse
import re

# DNA scoring matrix
MATRIX = ps.matrix_create("ACGT", 1, -1)
COMP = str.maketrans("ACGT", "TGCA")
VALID_QUERY = re.compile("[A|C|G|T]+")


def reverse_and_complement(query: str) -> str:
    return query.translate(COMP)[::-1]


def main() -> None:
    ag = argparse.ArgumentParser(
        prog="Pralign",
        description="Find a primer sequence that matches an arbitrarily long query DNA sequence",
    )
    ag.add_argument("query")
    ag.add_argument(
        "--reverse",
        "-r",
        help="Whether to reverse and make a complement of the query strand",
        action="store_true",
        default=False,
    )
    ag.add_argument(
        "--top", "-n", help="The top N primers to output into a csv file", default=20
    )
    ag.add_argument(
        "--longest",
        "-l",
        help="The length of the longest primer of interest",
        default=45,
    )
    ag.add_argument(
        "--shortest",
        "-s",
        help="The length of the shortest primer of interest",
        default=18,
    )
    namespace = ag.parse_args()
    query: str = namespace.query.strip().upper()
    if VALID_QUERY.fullmatch(query) is None:
        raise ArgumentError(
            f"query {query} is not valid. Only characters A, C, T, and G are allowed"
        )
    top: int = int(namespace.top)
    shortest: int = int(namespace.shortest)
    longest: int = int(namespace.longest)

    if namespace.reverse:
        query = reverse_and_complement(query)

    db: pathlib.Path = list(pathlib.Path(".").glob("*Primer*"))[
        -1
    ]  # select the last file in a lexicographically sorted list of primer databases.
    df: pl.DataFrame = (
        pl.read_excel(db)
        .lazy()
        .rename(
            {
                "Name *": "primer name",
                "Vial Number *": "vial number",
                "Sequence": "sequence",
            }
        )
        .select(["primer name", "vial number", "sequence"])
        .drop_nulls("sequence")  # some sequences are null
        .with_columns(pl.col("sequence").str.len_bytes().alias("sequence length"))
        .filter(
            (pl.col("sequence length") < longest)
            & (pl.col("sequence length") > shortest)
        )  # some sequences are question marks? Not relevant, but must set lower bound to at least 2
        .with_columns(
            pl.col("sequence")
            .map_elements(lambda seq: align(query, seq.upper()), return_dtype=pl.Int64)
            .alias("alignment score")
        )
        .collect()
    )
    df.sort(by="alignment score", descending=True).head(top).write_csv("alignments.csv")


def align(query: str, seq: str, gap_open: int = 1) -> int:
    result = ps.sw(query, seq, 1, 1, MATRIX)
    return result.score
