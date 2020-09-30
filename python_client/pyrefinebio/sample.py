import pyrefinebio.common.annotation as prb_annotation
import pyrefinebio.common.organism_index as prb_organism_index
import pyrefinebio.computational_result as prb_computational_result
import pyrefinebio.experiment as prb_experiment
import pyrefinebio.organism as prb_organism
import pyrefinebio.processor as prb_processor
from pyrefinebio.http import get_by_endpoint
from pyrefinebio.util import generator_from_pagination


class Sample:
    """Sample.

    get a sample based on accession code

        ex:
        >>> import pyrefinebio
        >>> accession_code = "GSM000000"
        >>> sample = pyrefinebio.Sample.get(accession_code)

    search for samples based on filters

        ex:
        >>> import pyrefinebio
        >>> samples = pyrefinebio.Sample.search(is_processed=True, specimen_part="soft-tissue sarcoma")
    """

    valid_filters = [
        "ordering",
        "title",
        "organism",
        "source_database",
        "source_archive_url",
        "has_raw",
        "platform_name",
        "technology",
        "manufacturer",
        "sex",
        "age",
        "specimen_part",
        "genotype",
        "disease",
        "disease_stage",
        "cell_line",
        "treatment",
        "race",
        "subject",
        "compound",
        "time",
        "is_processed",
        "is_public",
        "limit",
        "offset",
        "dataset_id",
        "experiment_accession_code",
        "accession_codes",
    ]

    def __init__(
        self=None,
        id=None,
        title=None,
        accession_code=None,
        source_database=None,
        organism=None,
        platform_accession_code=None,
        platform_name=None,
        pretty_platform=None,
        technology=None,
        manufacturer=None,
        protocol_info=None,
        annotations=[],
        results=[],
        source_archive_url=None,
        has_raw=None,
        sex=None,
        age=None,
        specimen_part=None,
        genotype=None,
        disease=None,
        disease_stage=None,
        cell_line=None,
        treatment=None,
        race=None,
        subject=None,
        compound=None,
        time=None,
        is_processed=None,
        created_at=None,
        last_modified=None,
        original_files=None,
        computed_files=None,
        experiment_accession_codes=None,
        experiments=None,
    ):
        self.id = id
        self.title = title
        self.accession_code = accession_code
        self.source_database = source_database
        self.organism = organism.Organism(**(organism or {}))
        self.platform_accession_code = platform_accession_code
        self.platform_name = platform_name
        self.pretty_platform = pretty_platform
        self.technology = technology
        self.manufacturer = manufacturer
        self.protocol_info = protocol_info
        self.annotations = [prb_annotation.Annotation(**annotation) for annotation in annotations]
        self.results = [
            prb_computational_result.ComputationalResult(**result) for result in results
        ]
        self.source_archive_url = source_archive_url
        self.has_raw = has_raw
        self.sex = sex
        self.age = age
        self.specimen_part = specimen_part
        self.genotype = genotype
        self.disease = disease
        self.disease_stage = disease_stage
        self.cell_line = cell_line
        self.treatment = treatment
        self.race = race
        self.subject = subject
        self.compound = compound
        self.time = time
        self.is_processed = is_processed
        self.created_at = created_at
        self.last_modified = last_modified
        self.original_files = original_files
        self.computed_files = computed_files
        self.experiment_accession_codes = experiment_accession_codes
        self.experiments = experiments

    @property
    def experiments(self):
        if not self._experiments:
            self._experiments = prb_experiment.Experiment.search(
                accession_code=self.experiment_accession_codes
            )

        return self._experiments

    @experiments.setter
    def experiments(self, value):
        self._experiments = value

    @classmethod
    def get(cls, accession_code):
        """Retrieve details about a sample based on its accession code.

        parameters:

            accession_code (str): The accession code for the sample to be retrieved.
        """

        result = get_by_endpoint("samples/" + accession_code)

        return cls(**result)

    @classmethod
    def search(cls, **kwargs):
        """Returns a list of samples.

        Parameters:
            ordering (str):                  which field to use when ordering the results
            title (str):
            organism (str):
            source_database (str):
            source_archive_url (str):
            has_raw (str):
            platform_name (str):
            technology (str):
            manufacturer (str):
            sex (str):
            age (number):
            specimen_part (str):
            genotype (str):
            disease (str):
            disease_stage (str):
            cell_line (str):
            treatment (str):
            race (str):
            subject (str):
            compound (str):
            time (str):
            is_processed (str):
            is_public (str):
            limit (int):                     Number of results to return per page.

            offset (int):                    The initial index from which to return the results.

            dataset_id (str):                Filters the result and only returns samples
                                             that are added to a dataset.

            experiment_accession_code (str): Filters the result and only returns only
                                             the samples associated with an experiment
                                             accession code.

            accession_codes (str):           Provide a list of sample accession codes
                                             separated by commas and the endpoint will
                                             only return information about these samples.
        """

        invalid_filters = []

        for filter in kwargs.keys():
            if filter not in cls.valid_filters:
                invalid_filters.append(filter)

        if invalid_filters:
            raise Exception("You supplied invalid filters - {0}".format(invalid_filters))

        result = get_by_endpoint("samples", params=kwargs)

        return generator_from_pagination(result, cls)
