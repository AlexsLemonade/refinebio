from bioinformatics_mill_models.models import Batch, SurveyJob
from .external_source import ExternalSourceSurveyor, ProcessorPipeline


class ArrayExpressSurveyor(ExternalSourceSurveyor):
    # Files API endpoint for ArrayExpress
    FILES_URL = "http://www.ebi.ac.uk/arrayexpress/json/v2/files"

    def source_type(self):
        return "ARRAY_EXPRESS"

    def determine_pipeline(self,
                           batch: Batch,
                           key_values: List[BatchKeyValue] = []):
        return ProcessorPipeline.MICRO_ARRAY_TO_PCL

    def survey(self, survey_job: SurveyJob):
        accession_code = 'A-AFFY-1'
        parameters = {'raw': 'true', 'array': platform}

        r = requests.get(FILES_URL, params=parameters)
        response_dictionary = r.json()

        try:
            experiments = response_dictionary['files']['experiment']
        except KeyError:  # If the platform does not exist or has no files...
            logger.info('No files were found with this platform accession code. ' +
                        'Try another accession code.')

        for experiment in experiments:
            data_files = experiment['file']

            # If there is only one file object in data_files, ArrayExpress does not
            # put it in a list of size 1 - This breaks the code if we attempt to
            # iterate over it like a list. The next section handles both cases.

            if (type(data_files) == list):
                for data_file in data_files:
                    if (data_file['kind'] == 'raw'):
                        url = data_file['url'].replace("\\", "")
                        self.handle_batch(Batch(size_in_bytes=0,
                                                download_url=url,
                                                raw_format="MICRO_ARRAY",
                                                accession_code=accession_code,
                                                organism=1))

            else:  # It is just one file object
                if (data_files['kind'] == 'raw'):
                    url = data_files['url'].replace("\\", "")
                    self.handle_batch(Batch(size_in_bytes=0,
                                            download_url=url,
                                            raw_format="MICRO_ARRAY",
                                            accession_code=accession_code,
                                            organism=1))

            break
