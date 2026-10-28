import requests as rq
from typing import Any

# NEB exposes an API
BASE_URI: str = "https://tmapi.neb.com"
# REQUEST_TYPE: str = (
#     "POST"  # The data is located in POST, but requests does not need it as an argument.
# )
# CONTENT_TYPE: str = "application/json"  # Content type must be application/json, but requests does it implicitly


def unnest(response: dict[str, Any]) -> Any:
    """Unnests a POST response.
    :param response: A JSON-like dict.
    """
    if isinstance(response, dict) and "data" in response.keys():
        response = response["data"]
    return response


# prodcodes:
# onetaq 2x std buf: onetaq-2
# phusion with HF buf: phusion-0
def get_tm_batch_onetaq_stdbuf(seqs: list[str]) -> Any:
    """Gets the Tm's for a batch of sequences for a PCR with OneTaq 2x MM with Standard Buffer.
    :param seqs: A list of sequences.
    :returns: Presumably a list of dictionaries containing information about each sequence, but should be verified.
    """
    r = rq.post(
        f"{BASE_URI}/tm/batch/",
        json={
            "seqpairs": [
                [seq] for seq in seqs
            ],  # seqpairs is a nested list, because it could take [seq1, seq2]
            "conc": 0.2,  # OneTaq protocol suggests 0.2 uM conc
            "prodcode": "onetaq-2",
        },
    )
    return unnest(r.json())


def get_tm_batch_phusion_hfbuf(seqs: list[str]) -> Any:
    """Gets the Tm's for a batch of sequences for a PCR with Phusion 2x MM with HF Buffer.
    :param seqs: A list of sequences.
    :returns: Presumably a list of dictionaries containing information about each sequence, but should be verified.
    """
    r = rq.post(
        f"{BASE_URI}/tm/batch/",
        json={
            "seqpairs": [[seq] for seq in seqs],
            "conc": 0.5,  # Phusion protocol suggests 0.5 uM conc
            "prodcode": "phusion-0",
        },
    )
    return unnest(r.json())


def __get_NEB_prodcodes():
    """A private function to retrieve the NEB product codes"""
    r = rq.get(f"{BASE_URI}/docs/productcodes")
    return r.text
