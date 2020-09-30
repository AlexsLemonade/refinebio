from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class Organism:
    def __init__(self, name=None, taxonomy_id=None):
        self.name = name
        self.taxonomy_id = taxonomy_id

    @classmethod
    def get(cls, name):
        response = get_by_endpoint("organisms/" + name)
        return Organism(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("organisms", params=kwargs)
        return generator_from_pagination(response, cls)
