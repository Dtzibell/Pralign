from typing import Any


def unnest(response: dict[str, Any]) -> Any:
    """Unnests a POST response.
    :param response: A JSON-like dict.
    """
    if isinstance(response, dict) and "data" in response.keys():
        response = response["data"]
    return response
