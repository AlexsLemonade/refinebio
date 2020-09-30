import pyrefinebio.computed_file as prb_computed_file
from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class Compendia:
    def __init__(
        self,
        id=None,
        primary_organism_name=None,
        organism_names=None,
        svd_algorithm=None,
        quant_sf_only=None,
        compendium_version=None,
        computed_file=None,
    ):
        self.id = id
        self.primary_organism_name = primary_organism_name
        self.organism_names = organism_names
        self.svd_algorithm = svd_algorithm
        self.quant_sf_only = quant_sf_only
        self.compendium_version = compendium_version
        self.computed_file = prb_computed_file.ComputedFile(**(computed_file or {}))

    @classmethod
    def get(cls, id):
        result = get_by_endpoint("compendia/" + str(id))
        return Compendia(**result)

    @classmethod
    def search(cls, **kwargs):
        result = get_by_endpoint("compendia", params=kwargs)
        return generator_from_pagination(result, cls)
