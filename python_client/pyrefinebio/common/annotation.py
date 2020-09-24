class Annotation:
    def __init__(self, **kwargs):
        self.data = kwargs["data"]
        self.is_ccdl = kwargs["is_ccdl"]
        self.created_at = kwargs["created_at"]
        self.last_modified = kwargs["last_modified"]
