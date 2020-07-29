import collections

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def requests_retry_session(
    retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504), session=None
):
    """
    Exponential back off for requests.

    via https://www.peterbe.com/plog/best-practice-with-retries-with-requests
    """
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


def flatten(d, parent_key="", sep="_"):
    """
    Flattens a dictionary using a seperator.

    via https://stackoverflow.com/a/6027615
    """
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def get_title_and_authors_for_pubmed_id(pmid):
    """ Given a PMID, return that PMID's (title, [authors]). """

    try:
        j_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id="
            + str(pmid)
            + "&retmode=json&tool=refinebio&email=hello@refine.bio"
        )
        resp = requests_retry_session().get(j_url, timeout=60)
        title = resp.json()["result"][str(pmid)]["title"]
        author_names = []
        for author in resp.json()["result"][str(pmid)]["authors"]:
            author_names.append(author["name"])

        return (title, author_names)
    except Exception:
        # This is fine for a timeout
        return ("", [])
