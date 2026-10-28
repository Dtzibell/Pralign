from ctypes import ArgumentError
import polars as pl
import parasail as ps
import pathlib
import argparse
import re
from .api import get_tm_batch_onetaq_stdbuf, get_tm_batch_phusion_hfbuf

# DNA scoring matrix
MATRIX = ps.matrix_create("ACGT", 1, -1)
COMP = str.maketrans("ACGT", "TGCA")
VALID_QUERY = re.compile("^[A|C|G|T]+$")


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
    ag.add_argument(
        "--output",
        "-o",
        help='The name of the output file, without the extension. E.g. "alignments"',
        default="alignments",
    )
    namespace = ag.parse_args()
    query: str = namespace.query.strip().upper()
    if VALID_QUERY.fullmatch(query) is None:
        raise ArgumentError(
            f"query {query} is not valid. Only characters A, C, T, and G are allowed"
        )
    top: int = int(namespace.top)
    shortest: int = int(namespace.shortest)
    if shortest < 2:
        raise ArgumentError(
            f'The argument "shortest" must be at least 2. Given {shortest}.'
        )
    longest: int = int(namespace.longest)
    output: str = namespace.output.strip()

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
