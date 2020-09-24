from pyrefinebio.common.annotation import Annotation
from pyrefinebio.common.organism_index import OgranismIndex
from pyrefinebio.http import get
from pyrefinebio.organism import Organism
from pyrefinebio.processor import Processor
from pyrefinebio.util import generator_from_pagination


class Sample:
    """
        Samples.

        get a sample based on accession code

            ex:
            >>> import pyrefinebio
            >>> accession_code = ""
            >>> sample = pyrefinebio.Samples.get(accession_code)

        search for sampes based on filters

            ex:
            >>> import pyrefinebio
            >>> samples = pyrefinebio.Samples.search(is_processed=True, specimen_part="soft-tissue sarcoma")
    """

    def __init__(self, **kwargs):
        self.id = kwargs["id"]
        self.title = kwargs["title"]
        self.accession_code = kwargs["accession_code"]
        self.source_database = kwargs["source_database"]
        self.organism = Organism(**kwargs["organism"])
        self.platform_accession_code = kwargs["platform_accession_code"]
        self.platform_name = kwargs["platform_name"]
        self.pretty_platform = kwargs["pretty_platform"]
        self.technology = kwargs["technology"]
        self.manufacturer = kwargs["manufacturer"]
        self.protocol_info = kwargs["protocol_info"]
        self.annotations = ([Annotation(**annotation) for annotation in kwargs["annotations"]],)
        self.results = [Result(**result) for result in kwargs["results"]]
        self.source_archive_url = kwargs["source_archive_url"]
        self.has_raw = kwargs["has_raw"]
        self.sex = kwargs["sex"]
        self.age = kwargs["age"]
        self.specimen_part = kwargs["specimen_part"]
        self.genotype = kwargs["genotype"]
        self.disease = kwargs["disease"]
        self.disease_stage = kwargs["disease_stage"]
        self.cell_line = kwargs["cell_line"]
        self.treatment = kwargs["treatment"]
        self.race = kwargs["race"]
        self.subject = kwargs["subject"]
        self.compound = kwargs["compound"]
        self.time = kwargs["time"]
        self.is_processed = kwargs["is_processed"]
        self.created_at = kwargs["created_at"]
        self.last_modified = kwargs["last_modified"]
        self.original_files = kwargs["original_files"]
        self.computed_files = kwargs["computed_files"]
        self.experiment_accession_codes = kwargs["computed_files"]

    @classmethod
    def get(cls, accession_code):
        """
            Retrieve details about a sample based on its accession code.

            parameters:

                accession_code (str): The accession code for the sample to be retrieved.
        """

        return cls(**get("samples/" + accession_code))

    @classmethod
    def search(
        cls,
        ordering="",
        title="",
        organism="",
        source_database="",
        source_archive_url="",
        has_raw="",
        platform_name="",
        technology="",
        manufacturer="",
        sex="",
        age=None,
        specimen_part="",
        genotype="",
        disease="",
        disease_stage="",
        cell_line="",
        treatment="",
        race="",
        subject="",
        compound="",
        time="",
        is_processed="",
        is_public="",
        limit=None,
        offset=None,
        dataset_id="",
        experiment_accession_code="",
        accession_codes="",
    ):
        """
        Returns a list of samples.

        Usage:
            `import pyrefinebio`
            `samples = pyrefinebio.Samples.llist()`

        Parameters:
            ordering (str):  which field to use when ordering the results
            title (str)
            organism (str)
            source_database (str)
            source_archive_url (str)
            has_raw (str)
            platform_name (str)
            technology (str)
            manufacturer (str)
            sex (str)
            age (number)
            specimen_part (str)
            genotype (str)
            disease (str)
            disease_stage (str)
            cell_line (str)
            treatment (str)
            race (str)
            subject (str)
            compound (str)
            time (str)
            is_processed (str)
            is_public (str)
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

        params = {
            "ordering": ordering,
            "title": title,
            "organism": organism,
            "source_database": source_database,
            "source_archive_url": source_archive_url,
            "has_raw": has_raw,
            "platform_name": platform_name,
            "technology": technology,
            "manufacturer": manufacturer,
            "sex": sex,
            "age": age,
            "specimen_part": specimen_part,
            "genotype": genotype,
            "disease": disease,
            "disease_stage": disease_stage,
            "cell_line": cell_line,
            "treatment": treatment,
            "race": race,
            "subject": subject,
            "compound": compound,
            "time": time,
            "is_processed": is_processed,
            "is_public": is_public,
            "limit": limit,
            "offset": offset,
            "dataset_id": dataset_id,
            "experiment_accession_code": experiment_accession_code,
            "accession_codes": accession_codes,
        }

        response = get("samples", params=params)

        return generator_from_pagination(response, cls)


class Result:
    def __init__(self, **kwargs):
        self.id = kwargs["id"]
        self.processor = Processor(**kwargs["processor"])

        if kwargs["organism_index"]:
            self.organism_index = OgranismIndex(**kwargs["organism_index"])
        else:
            self.organism_index = None
