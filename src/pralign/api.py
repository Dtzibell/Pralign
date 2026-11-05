import requests as rq
from typing import Any
from requests import ConnectionError

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


def batch_request(seqs: list[str], concentration: float, prodcode: str) -> Any:
    try:
        return unnest(
            rq.post(
                f"{BASE_URI}/tm/batch/",
                json={
                    "seqpairs": [
                        [seq] for seq in seqs
                    ],  # seqpairs is a nested list, because it could take [seq1, seq2]
                    "conc": 0.2,  # OneTaq protocol suggests 0.2 uM conc
                    "prodcode": "onetaq-2",
                },
                timeout=1,
            ).json()
        )
    except ConnectionError as e:
        # print(f"Caught exception {type(e).__name__}: {e}")
        lst = [{"tm1": None} for _ in seqs]
        return lst


# prodcodes:
# onetaq 2x std buf: onetaq-2
# phusion with HF buf: phusion-0
def get_tm_batch_onetaq_stdbuf(seqs: list[str]) -> Any:
    """Gets the Tm's for a batch of sequences for a PCR with OneTaq 2x MM with Standard Buffer.
    :param seqs: A list of sequences.
    :returns: Presumably a list of dictionaries containing information about each sequence, but should be verified.
    """
    prodcode = "onetaq-2"
    concentration = 0.2
    return batch_request(seqs, concentration, prodcode)
    # TODO: timeout


def get_tm_batch_phusion_hfbuf(seqs: list[str]) -> Any:
    """Gets the Tm's for a batch of sequences for a PCR with Phusion 2x MM with HF Buffer.
    :param seqs: A list of sequences.
    :returns: Presumably a list of dictionaries containing information about each sequence, but should be verified.
    """
    prodcode = "phusion-0"
    concentration = 0.5
    return batch_request(seqs, concentration, prodcode)


def __get_NEB_prodcodes():
    """A private function to retrieve the NEB product codes"""
    r = rq.get(f"{BASE_URI}/docs/productcodes")
    return r.text
