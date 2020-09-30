import pyrefinebio.computed_file as prb_computed_file
from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class Compendia:
    """Compendia.

    get a compendium result based on id

        ex:
        >>> import pyrefinebio
        >>> id = 1
        >>> sample = pyrefinebio.Compendia.get(id)

    search for compendium results based on filters

        ex:
        >>> import pyrefinebio
        >>> samples = pyrefinebio.Compendia.search(compendium_version="")
    """

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
        """Get a specific compendium result based on id

        parameters:

            id (int): the id for the compendium result you want to get
        """

        result = get_by_endpoint("compendia/" + str(id))
        return Compendia(**result)

    @classmethod
    def search(cls, **kwargs):
        """Search for a compendium result based on filters

        parameters:

            primary_organism__name (str):

            compendium_version (int):

            quant_sf_only (bool): True for RNA-seq Sample Compendium
                                  results or False for quantile normalized.

            result__id (int):

            ordering (str): Which field to use when ordering the results.

            limit (int): Number of results to return per page.

            offset (int): The initial index from which to return the results.

            latest_version (bool): True will only return the highest
                                   compendium_version for each primary_organism.
        """

        result = get_by_endpoint("compendia", params=kwargs)
        return generator_from_pagination(result, cls)
