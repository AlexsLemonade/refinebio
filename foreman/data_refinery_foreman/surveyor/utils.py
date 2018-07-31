import collections
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


def requests_retry_session(
    retries=3,
    backoff_factor=0.3,
    status_forcelist=(500, 502, 504),
    session=None,
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
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session

def flatten(d, parent_key='', sep='_'):
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

def get_title_for_pubmed_id(pmid):
    """ Given a PMID, return that PMID's title. """

    title = ""
    try:
        resp = requests.get("http://www.ncbi.nlm.nih.gov/pubmed/" + str(pmid), timeout=20)
        title = resp.text.split('<h1>')[1].split('</h1>')[0].strip()
    except:
        # This is fine for a timeout
        return ""

    return title
