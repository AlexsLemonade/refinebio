import requests
from typing import List, Dict

from data_refinery_common.models import (
    Batch,
    File,
    SurveyJobKeyValue,
    Organism,
)
from data_refinery_common.models.new_models import Experiment
from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


EXPERIMENTS_URL = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/"
SAMPLES_URL = EXPERIMENTS_URL + "{}/samples"
UNKNOWN="UNKNOWN"

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

    def get_experiment_metadata(self, experiment_accession_code: str) -> Dict:
        """

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

        # If there is more than one arraydesign listed in the experiment
        # then there is no other way to determine which array was used
        # for which sample other than looking at the header of the CEL
        # file. That obviously cannot happen until the CEL file has been
        # downloaded so we can just mark it as UNKNOWN and let the
        # downloader inspect the downloaded file to determine the
        # array then.
        if len(parsed_json["arraydesign"]) == 0:
            logger.warn("Experiment %s has no arraydesign listed.",
                        experiment_accession_code,
                        survey_job=self.survey_job.id)
            experiment["platform_accession_code"] = UNKNOWN
            experiment["platform_accession_name"] = UNKNOWN
        elif len(parsed_json["arraydesign"]) > 1:
            experiment["platform_accession_code"] = UNKNOWN
            experiment["platform_accession_name"] = UNKNOWN
        else:
            experiment["platform_accession_code"] = parsed_json["arraydesign"][0]["accession"]
            experiment["platform_accession_name"] = parsed_json["arraydesign"][0]["name"]

        experiment["release_date"] = parsed_json["releasedate"]

        if "lastupdatedate" in parsed_json:
            experiment["last_update_date"] = parsed_json["lastupdatedate"]
        else:
            experiment["last_update_date"] = parsed_json["releasedate"]

        try:
            experiment_object = Experiment.objects.get(accession_code=experiment_accession_code)
        except Experiment.DoesNotExist:
            experiment_object = Experiment()
            experiment_object.source_url = request_url
            experiment_object.name = parsed_json["name"]
            experiment_object.description = parsed_json["description"][0]["text"]
            experiment_object.platform_name = experiment["platform_accession_name"]
            experiment_object.platform_accession_code = experiment["platform_accession_code"]
            experiment_object.save()

        return experiment

    def _generate_samples(self,
                          samples: List[Dict],
                          experiment: Experiment
                          ) -> List[Sample]:
        """Generates a Sample item for each sample in an AE experiment..

        There are three possible data situations for a sample:

            - If the sample only has raw data available:
                Download this raw data and process it
            - If the sample has both raw and derived data:
                Download the raw data and process it
            - If the sample only has derived data:
                Download the derived data and no-op it.

        """
        for sample in samples:
            if "file" not in sample:
                continue

            organism_name = UNKNOWN
            for characteristic in sample["characteristic"]:
                if characteristic["category"].upper() == "ORGANISM":
                    organism_name = characteristic["value"].upper()

            if organism_name == UNKNOWN:
                logger.error("Sample from experiment %s did not specify the organism name.",
                             experiment["experiment_accession_code"],
                             survey_job=self.survey_job.id)
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

                file = File(name=sample_file["name"],
                            download_url=download_url,
                            raw_format=raw_format,
                            processed_format=processed_format,
                            size_in_bytes=-1)  # Will have to be determined later

                print ("Going to add batch..")
                self.add_batch(platform_accession_code=experiment["platform_accession_code"],
                               experiment_accession_code=experiment["experiment_accession_code"],
                               organism_id=organism_id,
                               organism_name=organism_name,
                               experiment_title=experiment["name"],
                               release_date=experiment["release_date"],
                               last_uploaded_date=experiment["last_update_date"],
                               files=[file])

    def discover_batches(self):

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

        experiment = self.get_experiment_metadata(experiment_accession_code)
        r = requests.get(SAMPLES_URL.format(experiment_accession_code))
        samples = r.json()["experiment"]["sample"]
        self._generate_samples(samples, experiment)
