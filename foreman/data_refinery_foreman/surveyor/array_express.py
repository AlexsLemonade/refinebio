import requests
from typing import List

from data_refinery_models.models import (
    Batch,
    BatchKeyValue,
    SurveyJob,
    SurveyJobKeyValue,
    Organism
)
from data_refinery_foreman.surveyor.external_source import (
    ExternalSourceSurveyor,
    ProcessorPipeline
)

# Import and set logger
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

EXPERIMENTS_URL = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/"
SAMPLES_URL = EXPERIMENTS_URL + "{}/samples"


class ArrayExpressSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return "ARRAY_EXPRESS"

    def downloader_task(self):
        return "data_refinery_workers.downloaders.array_express.download_array_express"

    def determine_pipeline(self,
                           batch: Batch,
                           key_values: List[BatchKeyValue] = []):
        return ProcessorPipeline.AFFY_TO_PCL

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

    def survey(self):
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

        for sample in samples:
            if "file" not in sample:
                continue

            organism_name = "UNKNOWN"
            for characteristic in sample["characteristic"]:
                if characteristic["category"].upper() == "ORGANISM":
                    organism_name = characteristic["value"].upper()

            if organism_name == "UNKNOWN":
                logger.error("Sample from experiment %s did not specify the organism name.",
                             experiment_accession_code)
                organism_id = 0
            else:
                organism_id = Organism.get_id_for_name(organism_name)

            for sample_file in sample["file"]:
                if sample_file["type"] != "data" and sample_file["name"] is not None:
                    continue

                self.handle_batch(Batch(
                    size_in_bytes=-1,  # Will have to be determined later
                    download_url=sample_file["url"],
                    raw_format=sample_file["name"].split(".")[-1],
                    processed_format="PCL",
                    platform_accession_code=experiment["platform_accession_code"],
                    experiment_accession_code=experiment_accession_code,
                    organism_id=organism_id,
                    organism_name=organism_name,
                    experiment_title=experiment["name"],
                    release_date=experiment["release_date"],
                    last_uploaded_date=experiment["last_update_date"],
                    name=sample_file["name"]
                ))
