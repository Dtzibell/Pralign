from ctypes import ArgumentError
import polars as pl
import parasail as ps
import pathlib
import argparse
import re
import requests as rq

# DNA scoring matrix
MATRIX = ps.matrix_create("ACGT", 1, -1)
COMP = str.maketrans("ACGT", "TGCA")
VALID_QUERY = re.compile("^[A|C|G|T]+$")
# NEB exposes an API
BASE_URI: str = "https://tmapi.neb.com"
TYPE: str = "POST"
CONTENT_TYPE: str = "application/json"


def unnest(response: dict):
    if isinstance(response, dict) and "data" in response.keys():
        response = response["data"]
    return response


# prodcodes:
# onetaq 2x std buf: onetaq-2
# phusion with HF buf: phusion-0
def get_batch_tm_onetaq_stdbuf(seqs: list[str]):
    r = rq.post(
        f"{BASE_URI}/tm/batch/",
        json={
            "seqpairs": [
                [seq] for seq in seqs
            ],  # seqpairs is a nested list, because it could take [seq1, seq2]
            "conc": 0.5,
            "prodcode": "onetaq-2",
        },
    )
    return unnest(r.json())


def get_batch_tm_phusion_hfbuf(seqs: list[str]):
    r = rq.post(
        f"{BASE_URI}/tm/batch/",
        json={
            "seqpairs": [[seq] for seq in seqs],
            "conc": 0.5,
            "prodcode": "phusion-0",
        },
    )
    return unnest(r.json())


def get_NEB_prodcodes():
    r = rq.get(f"{BASE_URI}/docs/productcodes")
    return r.text


def reverse_and_complement(query: str) -> str:
    return query.translate(COMP)[::-1]


def align(query: str, seq: str, gap_open: int = 1) -> dict[str, int | str]:
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
    longest: int = int(namespace.longest)
    output: str = namespace.output

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
        )
        .unnest("alignment")
        .filter(
            pl.col("aligned portion")
            .str.contains("-")
            .not_()
            .and_(pl.col("alignment score") > 16)
        )
        .sort(by="alignment score", descending=True)
        .with_columns(
            (pl.col("alignment score") / pl.col("sequence length")).alias(
                "aligned proportion"
            )
        )
        .sort(by="aligned proportion", descending=True)
        .head(top)
        .collect()
    )
    onetaq_tm = pl.Series(
        "OneTaq Tm",
        [
            dct["tm1"]
            for dct in get_batch_tm_onetaq_stdbuf(df["aligned portion"].to_list())
        ],
    )
    phusion_tm = pl.Series(
        "Phusion Tm",
        [
            dct["tm1"]
            for dct in get_batch_tm_phusion_hfbuf(df["aligned portion"].to_list())
        ],
    )
    df = df.with_columns(onetaq_tm, phusion_tm).drop("aligned portion")
    df.write_csv(f"{output}.csv", float_precision=3)
