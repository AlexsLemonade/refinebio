import pyrefinebio.organism as prb_organism
from pyrefinebio.http import get_by_endpoint


class QNTarget:
    def __init__(
        self,
        id=None,
        filename=None,
        size_in_bytes=None,
        is_qn_target=None,
        sha1=None,
        s3_bucket=None,
        s3_key=None,
        s3_url=None,
        created_at=None,
        last_modified=None,
        result=None,
    ):
        self.id = id
        self.filename = filename
        self.size_in_bytes = size_in_bytes
        self.is_qn_target = is_qn_target
        self.sha1 = sha1
        self.s3_bucket = s3_bucket
        self.s3_key = s3_key
        self.s3_url = s3_url
        self.created_at = created_at
        self.last_modified = last_modified
        self.result = result

    @classmethod
    def get(cls, organism_name):
        response = get_by_endpoint("qn_targets/" + organism_name)
        return QNTarget(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("qn_targets", params=kwargs)
        return [prb_organism.Organism(**qn_organism) for qn_organism in response]
