import logging

import requests

logger = logging.getLogger(__name__)


def get(endpoint, params=None):
    try:
        response = requests.get("http://localhost:8000/v1/" + endpoint + "/", params=params)

        response.raise_for_status()
        results = response.json()
        # cache results?

    except requests.exceptions.HTTPError:
        logging.error("GET %s failed", endpoint)
        return None

    return results


def post(endpoint, payload=None):
    try:
        response = requests.post(
            "http://localhost:8000/v1/" + endpoint + "/",
            data=payload,
            headers={"Content-Type": "application/json"},
        )

        response.raise_for_status()
        results = response.json()
        # cache results?

    except requests.exceptions.HTTPError:
        logging.error("POST %s failed", endpoint)
        return None

    return results


def put(endpoint, params=None, payload=None):
    try:
        response = requests.put(
            "http://localhost:8000/v1/" + endpoint + "/",
            data=payload,
            headers={"Content-Type": "application/json"},
        )

        response.raise_for_status()
        results = response.json()
        # cache results?

    except requests.exceptions.HTTPError:
        logging.error("POST %s failed", endpoint)
        return None

    return results
