class Processor:
    def __init__(self, id=None, name=None, version=None, docker_image=None, environment=None):
        self.id = id
        self.name = name
        self.version = version
        self.docker_image = docker_image
        self.environment = environment
