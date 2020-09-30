from pyrefinebio.http import get_by_endpoint


class Platform:
    def __init__(self, platform_accession_code=None, platform_name=None):
        self.platform_accession_code = platform_accession_code
        self.platform_name = platform_name

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("platforms", params=kwargs)
        return [Platform(**platform) for platform in response]
