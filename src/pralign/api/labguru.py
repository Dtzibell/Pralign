import os
import polars as pl
import requests as rq
from ctypes import ArgumentError

BASE_URI = "https://mftuebingen.labguru.com/api/v1/primers"


# loads the API_KEY
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
            f"{BASE_URI}",
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
