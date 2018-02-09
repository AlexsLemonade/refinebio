import requests
from typing import List, Dict

from data_refinery_common.models import (
    Batch,
    File,
    SurveyJobKeyValue,
    Organism,
)
from data_refinery_common.models.new_models import Experiment, Sample, ExperimentSampleAssociation
from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger

from django.conf import settings

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

        See an example at: https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-3050/samples
        Example:
        {
            'id': 511696,
            'accession': 'E-MTAB-3050',
            'name': 'Microarray analysis of in vitro differentiation of adult human pancreatic progenitor cells',
            'releasedate': '2014-10-31',
            'lastupdatedate': '2014-10-30',
            'organism': ['Homo sapiens'],
            'experimenttype': ['transcription profiling by array'],
            'experimentdesign': ['cell type comparison design',
                                 'development or differentiation design'],
            'description': [{'id': None,
                            'text': 'We have developed an in vitro culture system that allows us to expand progenitor cells from human islet preparations and differentiate them into insulin-producing cells. We noticed however, that cultures from individual islet preparations had very heterogeneous outcomes, from good differentiation to almost none. We therefore speculated that our progenitor cell cultures contained different kinds of cells and that the true endocrine progenitor cells are present in the successful cultures, but not in the unsuccessful ones. To address this issue and to begin to identify markers for the true endocrine progenitor cells we compared global gene expression between a very successful culture (final insulin-expression 10% of islets) and an unsuccessful one (final insulin-expression 0%). We also included RNA from freshly isolated islets for control purposes. The cultures from donor A yielded substantial differentiation, while the cultures from donor B showed no successful differentiation.'
                            }],
            'provider': [{'contact': 'Joel Habener', 'role': 'submitter',
                         'email': 'jhabener@partners.org'}],
            'samplecharacteristic': [{'category': 'age', 'value': ['38 year',
                                     '54 year']},
                                     {'category': 'developmental stage',
                                     'value': ['adult']},
                                     {'category': 'organism',
                                     'value': ['Homo sapiens']},
                                     {'category': 'organism part',
                                     'value': ['islet']}, {'category': 'sex',
                                     'value': ['female', 'male']}],
            'experimentalvariable': [{'name': 'cell type',
                                     'value': ['differentiated', 'expanded',
                                     'freshly isolated']}, {'name': 'individual'
                                     , 'value': ['A', 'B']},
                                     {'name': 'test result', 'value': ['NA',
                                     'successful', 'unsuccessful']}],
            'arraydesign': [{
                'id': 11048,
                'accession': 'A-AFFY-1',
                'name': 'Affymetrix GeneChip Human Genome U95Av2 [HG_U95Av2]',
                'count': 5,
                'legacy_id': 5728564,
                }],
            'protocol': [
                {'id': 1092859, 'accession': 'P-MTAB-41859'},
                {'id': 1092858, 'accession': 'P-MTAB-41860'},
                {'id': 1092857, 'accession': 'P-MTAB-41861'},
                {'id': 235758, 'accession': '  '},
                {'id': 1092863, 'accession': 'P-MTAB-41854'},
                {'id': 1092856, 'accession': 'P-MTAB-41862'},
                {'id': 1092860, 'accession': 'P-MTAB-41856'},
                {'id': 1092861, 'accession': 'P-MTAB-41855'},
                {'id': 1092862, 'accession': 'P-MTAB-41858'},
                {'id': 1092864, 'accession': 'P-MTAB-41857'},
                ],
            'bioassaydatagroup': [{
                'id': None,
                'name': 'rawData',
                'bioassaydatacubes': 5,
                'arraydesignprovider': None,
                'dataformat': 'rawData',
                'bioassays': 5,
                'isderived': 0,
                }, {
                'id': None,
                'name': 'processedData',
                'bioassaydatacubes': 5,
                'arraydesignprovider': None,
                'dataformat': 'processedData',
                'bioassays': 5,
                'isderived': 0,
                }, {
                'id': None,
                'name': 'image',
                'bioassaydatacubes': 5,
                'arraydesignprovider': None,
                'dataformat': 'image',
                'bioassays': 5,
                'isderived': 0,
                }],
            }

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
            experiment_object.name = parsed_json["name"]
            experiment_object.description = parsed_json["description"][0]["text"]
            experiment_object.platform_name = experiment["platform_accession_name"]
            experiment_object.platform_accession_code = experiment["platform_accession_code"]
            experiment_object.save()

            # TODO: Is there other K/V pair data we should create here?

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

            sample_accession_code = sample["assay"]["name"]

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
                filename = None

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

                if not download_url:
                    logger.error("Sample %s from experiment %s did not specify a download url, skipping.",
                            sample_accession_code,
                            experiment.accession_code,
                            survey_job=self.survey_job.id)
                    continue
                filename = sub_file['name']

            if not download_url:
                print(sub_file)

            # Create the sample object
            try:
                sample = Sample.objects.get(accession_code=sample_accession_code)
                logger.error("Sample %s from experiment %s already exists, skipping object creation.",
                         sample_accession_code,
                         experiment.accession_code,
                         survey_job=self.survey_job.id)
                continue
            except Sample.DoesNotExist:
                sample = Sample()
                sample.accession_code = sample_accession_code
                sample.source_archive_url = download_url
                sample.source_filename = filename
                sample.is_downloaded = False
                sample.has_raw = has_raw
                sample.organism = organism
                sample.save()

            association = ExperimentSampleAssociation()
            association.experiment = experiment
            association.sample = sample
            association.save()

            logger.info("Created Sample: " + str(sample))

            created_samples.append(sample)

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
            return []

        samples = self.create_samples_from_api(experiment)
        return experiment, samples
