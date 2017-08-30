import requests
from typing import List, Dict

from data_refinery_models.models import (
    Batch,
    BatchKeyValue,
    SurveyJobKeyValue,
    Organism
)
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders


# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


EXPERIMENTS_URL = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/"
SAMPLES_URL = EXPERIMENTS_URL + "{}/samples"


class ArrayExpressSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return Downloaders.ARRAY_EXPRESS.value

    def determine_pipeline(self,
                           batch: Batch,
                           key_values: List[BatchKeyValue] = []):
        # If it's a CEL file run SCAN.UPC on it.
        if batch.raw_format == "CEL":
            return ProcessorPipeline.AFFY_TO_PCL
        # If only processed data is available then we don't need to do
        # anything to it
        elif batch.raw_format == batch.processed_format:
            return ProcessorPipeline.NO_OP
        # If it's not CEL and it's not already processed then we just
        # want to download it for Jackie's grant.
        else:
            return ProcessorPipeline.NONE

    @staticmethod
    def get_experiment_metadata(experiment_accession_code):
        experiment_request = requests.get(EXPERIMENTS_URL + experiment_accession_code)
        parsed_json = experiment_request.json()["experiments"]["experiment"][0]

        experiment = {}

        experiment["name"] = parsed_json["name"]
        experiment["experiment_accession_code"] = experiment_accession_code

        # If there is more than one arraydesign listed in the experiment
        # then there is no other way to determine which array was used
        # for which sample other than looking at the header of the CEL
        # file. That obviously cannot happen until the CEL file has been
        # downloaded so we can just mark it as UNKNOWN and let the
        # downloader inspect the downloaded file to determine the
        # array then.
        if len(parsed_json["arraydesign"]) == 0:
            logger.warn("Experiment %s has no arraydesign listed.", experiment_accession_code)
            experiment["platform_accession_code"] = "UNKNOWN"
        elif len(parsed_json["arraydesign"]) > 1:
            experiment["platform_accession_code"] = "UNKNOWN"
        else:
            experiment["platform_accession_code"] = \
                parsed_json["arraydesign"][0]["accession"]

        experiment["release_date"] = parsed_json["releasedate"]

        if "lastupdatedate" in parsed_json:
            experiment["last_update_date"] = parsed_json["lastupdatedate"]
        else:
            experiment["last_update_date"] = parsed_json["releasedate"]

        return experiment

    def _generate_batches(self,
                          samples: List[Dict],
                          experiment: Dict,
                          replicate_raw: bool = True) -> List[Batch]:
        """Generates a Batch for each sample in samples.

        Uses the metadata contained in experiment (which should be
        generated via get_experiment_metadata) to add additional
        metadata to each Batch. If replicate_raw is True (the default)
        then only raw files will be replicated. Otherwise all files
        will be replicated.
        """
        batches = []
        for sample in samples:
            if "file" not in sample:
                continue

            organism_name = "UNKNOWN"
            for characteristic in sample["characteristic"]:
                if characteristic["category"].upper() == "ORGANISM":
                    organism_name = characteristic["value"].upper()

            if organism_name == "UNKNOWN":
                logger.error("Sample from experiment %s did not specify the organism name.",
                             experiment["experiment_accession_code"])
                organism_id = 0
            else:
                organism_id = Organism.get_id_for_name(organism_name)

            for sample_file in sample["file"]:
                # Generally we only want to replicate the raw data if
                # we can, however if there isn't raw data then we'll
                # take the processed stuff.
                if (replicate_raw and sample_file["type"] != "data") \
                        or sample_file["name"] is None:
                    continue

                # sample_file["comment"] is only a list if there's
                # more than one comment...
                comments = sample_file["comment"]
                if isinstance(comments, list):
                    # Could be: "Derived ArrayExpress Data Matrix FTP
                    # file" or: "ArrayExpress FTP file". If there is
                    # no comment with a name including "FTP file" then
                    # we don't know where to download it so we need to
                    # mark this job as an error. Therefore don't catch
                    # the potential exception where download_url
                    # doesn't get defined.
                    for comment in comments:
                        if comment["name"].find("FTP file") != -1:
                            download_url = comment["value"]
                else:
                    download_url = comments["value"]

                raw_format = sample_file["name"].split(".")[-1]
                processed_format = "PCL" if replicate_raw else raw_format

                batches.append(Batch(
                    size_in_bytes=-1,  # Will have to be determined later
                    download_url=download_url,
                    raw_format=raw_format,
                    processed_format=processed_format,
                    platform_accession_code=experiment["platform_accession_code"],
                    experiment_accession_code=experiment["experiment_accession_code"],
                    organism_id=organism_id,
                    organism_name=organism_name,
                    experiment_title=experiment["name"],
                    release_date=experiment["release_date"],
                    last_uploaded_date=experiment["last_update_date"],
                    name=sample_file["name"]
                ))

        return batches

    def survey(self) -> bool:
        experiment_accession_code = (
            SurveyJobKeyValue
            .objects
            .get(survey_job_id=self.survey_job.id,
                 key__exact="experiment_accession_code")
            .value
        )

        logger.info("Surveying experiment with accession code: %s.", experiment_accession_code)

        experiment = self.get_experiment_metadata(experiment_accession_code)
        r = requests.get(SAMPLES_URL.format(experiment_accession_code))
        samples = r.json()["experiment"]["sample"]
        batches = self._generate_batches(samples, experiment)

        if len(samples) != 0 and len(batches) == 0:
            # Found no samples with raw data, so replicate the
            # processed data instead
            batches = self._generate_batches(samples, experiment, replicate_raw=False)

        # Group batches based on their download URL and handle each group.
        download_urls = {batch.download_url for batch in batches}
        for url in download_urls:
            batches_with_url = [batch for batch in batches if batch.download_url == url]
            try:
                self.handle_batches(batches_with_url)
            except Exception:
                return False

        return True
