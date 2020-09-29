class OrganismIndex:
    def __init__(
        self,
        id=None,
        assembly_name=None,
        organism_name=None,
        source_version=None,
        index_type=None,
        salmon_version=None,
        download_url=None,
        result_id=None,
        last_modified=None,
    ):
        self.id = id
        self.assembly_name = assembly_name
        self.organism_name = organism_name
        self.source_version = source_version
        self.index_type = index_type
        self.salmon_version = salmon_version
        self.download_url = download_url
        self.result_id = result_id
        self.last_modified = last_modified
