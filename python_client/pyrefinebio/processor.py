from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class Processor:
    def __init__(self, id=None, name=None, version=None, docker_image=None, environment=None):
        self.id = id
        self.name = name
        self.version = version
        self.docker_image = docker_image
        self.environment = environment

    @classmethod
    def get(cls, id):
        response = get_by_endpoint("processors/" + str(id))
        return Processor(**response)

    @classmethod
    def search(cls, **kwargs):
        response = get_by_endpoint("processors", params=kwargs)
        return generator_from_pagination(response, cls)
