import logging
import os

import requests

logger = logging.getLogger(__name__)

base_url = os.getenv("BASE_URL") or "https://api.refine.bio/v1/"


def get(url, params=None):
    response = requests.get(url, params=params)

    response.raise_for_status()

    return response.json()


def post(url, payload=None):
    response = requests.post(url, data=payload, headers={"Content-Type": "application/json"})

    response.raise_for_status()

    return response.json()


def put(url, payload=None):
    response = requests.put(url, data=payload, headers={"Content-Type": "application/json"})

    response.raise_for_status()

    return response.json()


def get_by_endpoint(endpoint, params=None):
    url = base_url + endpoint + "/"

    return get(url, params=params)


def post_by_endpoint(endpoint, payload=None):
    url = base_url + endpoint + "/"

    return post(url, payload=payload)


def put_by_endpoint(endpoint, payload=None):
    url = base_url + endpoint + "/"

    return put(url, payload=payload)
