class OgranismIndex:
    def __init__(self, **kwargs):
        self.id = kwargs["id"]
        self.assembly_name = kwargs["assembly_name"]
        self.organism_name = kwargs["organism_name"]
        self.source_version = kwargs["source_version"]
        self.index_type = kwargs["index_type"]
        self.salmon_version = kwargs["salmon_version"]
        self.download_url = kwargs["download_url"]
        self.result_id = kwargs["result_id"]
        self.last_modified = kwargs["last_modified"]
