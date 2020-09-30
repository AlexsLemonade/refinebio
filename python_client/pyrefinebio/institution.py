from pyrefinebio.http import get_by_endpoint


class Institution:
    def __init__(self, submitter_institution=None):
        self.submitter_institution = submitter_institution

    @classmethod
    def search(cls):
        response = get_by_endpoint("institutions")
        return [Institution(**institution) for institution in response]
