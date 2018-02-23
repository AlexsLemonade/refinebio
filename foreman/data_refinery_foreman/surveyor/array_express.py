import requests

from django.conf import settings
from typing import List, Dict

from data_refinery_common.models import (
    Batch,
    File,
    SurveyJobKeyValue,
    Organism,
    Experiment,
    ExperimentAnnotation,
    Sample,
    SampleAnnotation,
    ExperimentSampleAssociation,
    OriginalFile
)
from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)


EXPERIMENTS_URL = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/"
SAMPLES_URL = EXPERIMENTS_URL + "{}/samples"
UNKNOWN="UNKNOWN"

class UnsupportedPlatformException(BaseException):
    pass

class ArrayExpressSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return Downloaders.ARRAY_EXPRESS.value

    def determine_pipeline(self, batch: Batch, key_values: Dict = {}) -> ProcessorPipeline:
        # If it's a CEL file run SCAN.UPC on it.
        if batch.files[0].raw_format == "CEL":
            return ProcessorPipeline.AFFY_TO_PCL
        # If only processed data is available then we don't need to do
        # anything to it
        elif batch.files[0].raw_format == batch.files[0].processed_format:
            return ProcessorPipeline.NO_OP
        # If it's not CEL and it's not already processed then we just
        # want to download it for Jackie's grant.
        else:
            return ProcessorPipeline.NONE

    def group_batches(self) -> List[List[Batch]]:
        return utils.group_batches_by_first_file(self.batches)

    def create_experiment_from_api(self, experiment_accession_code: str) -> Dict:
        """
        Given an experiment accession code, create an Experiment object.

        Will raise an UnsupportedPlatformException if this experiment was
        conducting using a platform which we don't support.

        See an example at: https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-3050/sample
        """
        request_url = EXPERIMENTS_URL + experiment_accession_code;
        experiment_request = requests.get(request_url)
        parsed_json = experiment_request.json()["experiments"]["experiment"][0]

        experiment = {}

        experiment["name"] = parsed_json["name"]
        experiment["experiment_accession_code"] = experiment_accession_code

        # This experiment has no platform at all, and is therefore useless.
        if 'arraydesign' not in parsed_json or len(parsed_json["arraydesign"]) == 0:
            logger.warn("Experiment %s has no arraydesign listed.",
                        experiment_accession_code,
                        survey_job=self.survey_job.id)
            raise UnsupportedPlatformException
        # If there is more than one arraydesign listed in the experiment
        # then there is no other way to determine which array was used
        # for which sample other than looking at the header of the CEL
        # file. That obviously cannot happen until the CEL file has been
        # downloaded so we can just mark it as UNKNOWN and let the
        # downloader inspect the downloaded file to determine the
        # array then.
        elif len(parsed_json["arraydesign"]) > 1:
            experiment["platform_accession_code"] = UNKNOWN
            experiment["platform_accession_name"] = UNKNOWN
        else:
            if parsed_json["arraydesign"][0]["accession"] not in settings.SUPPORTED_PLATFORMS:
                raise UnsupportedPlatformException

            experiment["platform_accession_code"] = parsed_json["arraydesign"][0]["accession"]
            experiment["platform_accession_name"] = parsed_json["arraydesign"][0]["name"]

        experiment["release_date"] = parsed_json["releasedate"]

        if "lastupdatedate" in parsed_json:
            experiment["last_update_date"] = parsed_json["lastupdatedate"]
        else:
            experiment["last_update_date"] = parsed_json["releasedate"]

        # Create the experiment object
        try:
            experiment_object = Experiment.objects.get(accession_code=experiment_accession_code)
            logger.error("Experiment %s already exists, skipping object creation.",
                experiment_accession_code,
                survey_job=self.survey_job.id)
        except Experiment.DoesNotExist:
            experiment_object = Experiment()
            experiment_object.accession_code = experiment_accession_code
            experiment_object.source_url = request_url
            experiment_object.source_database = "ARRAY_EXPRESS"
            experiment_object.name = parsed_json["name"]
            experiment_object.description = parsed_json["description"][0]["text"]
            experiment_object.platform_name = experiment["platform_accession_name"]
            experiment_object.platform_accession_code = experiment["platform_accession_code"]
            experiment_object.save()

            json_xa = ExperimentAnnotation()
            json_xa.experiment = experiment_object
            json_xa.data = parsed_json
            json_xa.is_ccdl = False
            json_xa.save()

            ## Fetch and parse the IDF file for any other fields
            IDF_URL_TEMPLATE = "https://www.ebi.ac.uk/arrayexpress/files/{code}/{code}.idf.txt"
            idf_url = IDF_URL_TEMPLATE.format(code=experiment_accession_code)
            idf_text = requests.get(idf_url).text
            lines = idf_text.split('\n')
            idf_dict = {}
            for line in lines:
                keyval = line.strip().split('\t')
                if len(keyval) == 2:
                    idf_dict[keyval[0]] = keyval[1]
                elif len(keyval) > 2:
                    idf_dict[keyval[0]] = keyval[1:]

            idf_xa = ExperimentAnnotation()
            idf_xa.data = idf_dict
            idf_xa.experiment = experiment_object
            idf_xa.is_ccdl = False
            idf_xa.save()

        return experiment_object

    def create_samples_from_api(self,
                          experiment: Experiment
                          ) -> List[Sample]:
        """Generates a Sample item for each sample in an AE experiment..

        There are many possible data situations for a sample:

            - If the sample only has raw data available:
                - If it is on a platform that we support:
                    Download this raw data and process it
                - If it is not on a platform we support:
                    Don't download anything, don't process anything
            - If the sample has both raw and derived data:
                - If the raw data is on a platform we support:
                    Download the raw data and process it, abandon the derived data
                - If the raw data is not on a platform we support
                    Download the derived data and no-op it, abandon the raw data
            - If the sample only has derived data:
                Download the derived data and no-op it.

        See an example at: https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-3050/samples

        """

        created_samples = []

        r = requests.get(SAMPLES_URL.format(experiment.accession_code))
        samples = r.json()["experiment"]["sample"]

        # An experiment can have many samples
        for sample in samples:

            # For some reason, this sample has no files associated with it.
            if "file" not in sample:
                continue

            # XXX: Somebody needs to explain this to me.
            try:
                sample_accession_code = sample["assay"]["name"]
            except KeyError:
                sample_accession_code = sample["source"]["name"]

            # Figure out the Organism for this sample
            organism_name = UNKNOWN
            for characteristic in sample["characteristic"]:
                if characteristic["category"].upper() == "ORGANISM":
                    organism_name = characteristic["value"].upper()

            if organism_name == UNKNOWN:
                logger.error("Sample from experiment %s did not specify the organism name.",
                             experiment.accession_code,
                             survey_job=self.survey_job.id)
                organism = None
            else:
                organism = Organism.get_object_for_name(organism_name)

            # A sample may actually have many sub files.
            # If there is raw data, take that.
            # If not, take the derived.
            has_raw = False
            for sub_file in sample['file']:
                # Some have the 'data' field, but not the actual data
                # Ex: E-GEOD-9656
                if sub_file['type'] == "data" and sub_file['comment']['value'] != None:
                    has_raw = True

            for sub_file in sample['file']:

                # Skip derived data if we have it raw.
                if has_raw:
                    if sub_file['type'] == "derived data":
                        continue

                download_url = None
                filename = sub_file["name"]

                # sub_file["comment"] is only a list if there's
                # more than one comment...
                comments = sub_file["comment"]
                if isinstance(comments, list):
                    # Could be: "Derived ArrayExpress Data Matrix FTP
                    # file" or: "ArrayExpress FTP file". If there is
                    # no comment with a name including "FTP file" then
                    # we don't know where to download it so we need to
                    # mark this job as an error. Therefore don't catch
                    # the potential exception where download_url
                    # doesn't get defined.
                    for comment in comments:
                        if "FTP file" in comment["name"]:
                            download_url = comment["value"]
                            break
                else:
                    download_url = comments["value"]
                    filename = sub_file["name"]

                if not download_url:
                    logger.error("Sample %s from experiment %s did not specify a download url, skipping.",
                            sample_accession_code,
                            experiment.accession_code,
                            survey_job=self.survey_job.id)
                    continue

            if not download_url:
                print(sub_file)

            # Create the sample object
            try:
                sample_object = Sample.objects.get(accession_code=sample_accession_code)
                logger.error("Sample %s from experiment %s already exists, skipping object creation.",
                         sample_accession_code,
                         experiment.accession_code,
                         survey_job=self.survey_job.id)
                continue
            except Sample.DoesNotExist:
                sample_object = Sample()
                sample_object.accession_code = sample_accession_code
                sample_object.organism = organism
                sample_object.save()

                # This creates key values for a given sample.
                # This looks a bit of a mess.

                # XXX: TODO: Harmonize desired values to filter on and 
                # extact to properties of the Sample itself.

                sample_annotation = SampleAnnotation()
                sample_annotation.data = sample
                sample.sample = sample_object
                sample.is_ccdl = False
                sample.save()

                original_file = OriginalFile()
                original_file.sample = sample_object
                original_file.source_filename = filename
                original_file.source_url = download_url
                original_file.is_downloaded = False
                original_file.is_archive = True
                original_file.has_raw = has_raw
                original_file.save()

            association = ExperimentSampleAssociation()
            association.experiment = experiment
            association.sample = sample_object
            association.save()

            logger.info("Created Sample: " + str(sample_object))

            created_samples.append(sample_object)

        return created_samples

    def discover_experiment_and_samples(self) -> (Experiment, List[Sample]):

        experiment_accession_code = (
            SurveyJobKeyValue
            .objects
            .get(survey_job_id=self.survey_job.id,
                 key__exact="experiment_accession_code")
            .value
        )

        logger.info("Surveying experiment with accession code: %s.",
                    experiment_accession_code,
                    survey_job=self.survey_job.id)

        try:
            experiment = self.create_experiment_from_api(experiment_accession_code)
        except UnsupportedPlatformException as e:
            logger.info("Experiment with accession code: %s was not on a supported platform, skipping.",
                experiment_accession_code,
                survey_job=self.survey_job.id)
            return None, []

        samples = self.create_samples_from_api(experiment)
        return experiment, samples
