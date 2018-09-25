from typing import List
from data_refinery_common.models import Batch


def group_batches_by_first_file(batches: List[Batch]) -> List[List[Batch]]:
    """Groups batches based on the download URL of their first File."""

    # Builds a mapping of each unique download_url to a list of
    # Batches whose first File's download_url matches.
    download_url_mapping = {}
    for batch in batches:
        download_url = batch.files[0].download_url
        if download_url in download_url_mapping:
            download_url_mapping[download_url].append(batch)
        else:
            download_url_mapping[download_url] = [batch]

    # The values of the mapping we built are the groups the
    # batches should be grouped into.
    return list(download_url_mapping.values())
