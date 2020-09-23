from pyrefinebio.http import get


class Samples:
    """
        Samples.
    """

    def get(self, accession_code):
        """
            Retrieve details about a sample based on its accession code.

            parameters:

                accession_code (str): The accession code for the sample to be retrieved.
        """

        return get("samples/" + accession_code)

    def llist(
        self,
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
        age=0,
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
        limit=0,
        offset=0,
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
            ordering: "",
            title: "",
            organism: "",
            source_database: "",
            source_archive_url: "",
            has_raw: "",
            platform_name: "",
            technology: "",
            manufacturer: "",
            sex: "",
            age: 0,
            specimen_part: "",
            genotype: "",
            disease: "",
            disease_stage: "",
            cell_line: "",
            treatment: "",
            race: "",
            subject: "",
            compound: "",
            time: "",
            is_processed: "",
            is_public: "",
            limit: 0,
            offset: 0,
            dataset_id: "",
            experiment_accession_code: "",
            accession_codes: "",
        }

        return get("samples", params=params)
        # TODO: iterate through pagination
        # do some stuff with generators?
