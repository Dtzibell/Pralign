import polars as pl
import parasail as ps
import pathlib
import argparse

# DNA scoring matrix
MATRIX = ps.matrix_create("ACGT", 1, -1)
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_and_complement(query: str) -> str:
    return query.translate(COMP)[::-1]


def main() -> None:
    ag = argparse.ArgumentParser()
    ag.add_argument("query")
    ag.add_argument("--reverse", "-r", action="store_true")
    namespace = ag.parse_args()
    query = namespace.query.upper()
    if namespace.reverse:
        query = reverse_and_complement(query)

    db = list(pathlib.Path(".").glob("primers*"))[
        -1
    ]  # select the most recent primers file
    df = (
        pl.read_excel(db)
        .lazy()
        .select(["Name *", "Vial Number *", "Sequence"])
        .drop_nulls("Sequence")
        .with_columns(pl.col("Sequence").str.len_bytes().alias("seq_len"))
        .filter((pl.col("seq_len") < 40) & (pl.col("seq_len") > 18))
        .sort(by="seq_len", descending=False)
        .with_columns(
            pl.col("Sequence")
            .map_elements(lambda seq: align(query, seq.upper()), return_dtype=pl.Int64)
            .alias("align_score")
        )
        .collect()
    )
    print(df.height, df["seq_len"].min(), df["seq_len"].max())
    df = df.sort(by="align_score", descending=True)
    df.head(20).write_csv("alignments.csv")


def align(query: str, seq: str, gap_open: int = 1) -> int:
    result = ps.sw(query, seq, 1, 1, MATRIX)
    return result.score
