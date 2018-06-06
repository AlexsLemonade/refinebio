import requests

from django.utils.dateparse import parse_datetime
from typing import List, Dict

from data_refinery_common.models import (
    SurveyJobKeyValue,
    Organism,
    Experiment,
    ExperimentAnnotation,
    Sample,
    SampleAnnotation,
    ExperimentSampleAssociation,
    OriginalFile,
    OriginalFileSampleAssociation,
    ExperimentOrganismAssociation
)

from data_refinery_foreman.surveyor import harmony, utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.utils import (
    get_supported_microarray_platforms,
    get_readable_platform_names)

logger = get_and_configure_logger(__name__)


EXPERIMENTS_URL = "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/"
SAMPLES_URL = EXPERIMENTS_URL + "{}/samples"
UNKNOWN = "UNKNOWN"


class UnsupportedPlatformException(Exception):
    pass


class ArrayExpressSurveyor(ExternalSourceSurveyor):

    def source_type(self):
        return Downloaders.ARRAY_EXPRESS.value

    def create_experiment_from_api(self, experiment_accession_code: str) -> (Experiment, Dict):
        """Given an experiment accession code, create an Experiment object.

        Also returns a dictionary of additional information about the
        platform discovered for the experiment.

        Will raise an UnsupportedPlatformException if this experiment was
        conducting using a platform which we don't support.

        See an example at: https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/E-MTAB-3050/sample
        """
        request_url = EXPERIMENTS_URL + experiment_accession_code
        experiment_request = utils.requests_retry_session().get(request_url, timeout=60)

        try:
            parsed_json = experiment_request.json()["experiments"]["experiment"][0]
        except KeyError:
            logger.error("Remote experiment has no Experiment data!",
                         experiment_accession_code=experiment_accession_code,
                         survey_job=self.survey_job.id)
            raise

        experiment = {}
        platform_warning = False

        experiment["name"] = parsed_json["name"]
        experiment["experiment_accession_code"] = experiment_accession_code

        # This experiment has no platform at all, and is therefore useless.
        if 'arraydesign' not in parsed_json or len(parsed_json["arraydesign"]) == 0:
            logger.warn("Remote experiment has no arraydesign listed.",
                        experiment_accession_code=experiment_accession_code,
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
            external_accession = parsed_json["arraydesign"][0]["accession"]
            for platform in get_supported_microarray_platforms():
                if platform["external_accession"] == external_accession:
                    experiment["platform_accession_code"] = platform["platform_accession"]
                    platform_mapping = get_readable_platform_names()
                    experiment["platform_accession_name"] = platform_mapping[
                        platform["platform_accession"]]

            if "platform_accession_code" not in experiment:
                # We don't know what platform this accession corresponds to.
                experiment["platform_accession_code"] = external_accession
                experiment["platform_accession_name"] = UNKNOWN
                platform_warning = True

        experiment["release_date"] = parsed_json["releasedate"]

        if "lastupdatedate" in parsed_json:
            experiment["last_update_date"] = parsed_json["lastupdatedate"]
        else:
            experiment["last_update_date"] = parsed_json["releasedate"]

        # Create the experiment object
        try:
            experiment_object = Experiment.objects.get(accession_code=experiment_accession_code)
            logger.debug("Experiment already exists, skipping object creation.",
                         experiment_accession_code=experiment_accession_code,
                         survey_job=self.survey_job.id)
        except Experiment.DoesNotExist:

            experiment_object = Experiment()
            experiment_object.accession_code = experiment_accession_code
            experiment_object.source_url = request_url
            experiment_object.source_database = "ARRAY_EXPRESS"
            experiment_object.name = parsed_json["name"]
            experiment_object.technology = "MICROARRAY"
            experiment_object.description = parsed_json["description"][0]["text"]
            experiment_object.source_first_published = parse_datetime(experiment["release_date"])
            experiment_object.source_last_modified = parse_datetime(experiment["last_update_date"])
            experiment_object.save()

            # We still create the Experiment and samples if there is processed
            # data but we can't support the reprocessing of raw data.
            if platform_warning:
                warning_xa = ExperimentAnnotation()
                warning_xa.experiment = experiment_object
                warning_xa.data = {'has_unsupported_platform': True}
                warning_xa.is_ccdl = False
                warning_xa.save()

            json_xa = ExperimentAnnotation()
            json_xa.experiment = experiment_object
            json_xa.data = parsed_json
            json_xa.is_ccdl = False
            json_xa.save()

            ## Fetch and parse the IDF/SDRF file for any other fields
            IDF_URL_TEMPLATE = "https://www.ebi.ac.uk/arrayexpress/files/{code}/{code}.idf.txt"
            SDRF_URL_TEMPLATE = "https://www.ebi.ac.uk/arrayexpress/files/{code}/{code}.sdrf.txt"
            idf_url = IDF_URL_TEMPLATE.format(code=experiment_accession_code)
            sdrf_url = SDRF_URL_TEMPLATE.format(code=experiment_accession_code)
            idf_text = utils.requests_retry_session().get(idf_url, timeout=60).text

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

            if 'Investigation Title' in idf_dict:
                experiment_object.title = idf_dict['Investigation Title']
            if 'Person Affiliation' in idf_dict:
                # This is very rare, ex: E-MEXP-32
                if isinstance(idf_dict['Person Affiliation'], list):
                    unique_people = list(set(idf_dict['Person Affiliation']))
                    experiment_object.submitter_institution = ", ".join(unique_people)[:255]
                else:
                    experiment_object.submitter_institution = idf_dict['Person Affiliation']
            if 'Protocol Description' in idf_dict:
                experiment_object.protocol_description = ", ".join(idf_dict['Protocol Description'])
            if 'Publication Title' in idf_dict:
                # This will happen for some superseries.
                # Ex: E-GEOD-29536
                # Assume most recent is "best:, store the rest in experiment annotation.
                if isinstance(idf_dict['Publication Title'], list):
                    experiment_object.publication_title = idf_dict['Publication Title'][0]
                else:
                    experiment_object.publication_title = idf_dict['Publication Title']
                experiment_object.has_publication = True
            if 'Publication DOI' in idf_dict:
                if isinstance(idf_dict['Publication DOI'], list):
                    experiment_object.publication_doi = idf_dict['Publication DOI'][0]
                else:
                    experiment_object.publication_doi = idf_dict['Publication DOI']
                experiment_object.has_publication = True
            if 'PubMed ID' in idf_dict:
                if isinstance(idf_dict['PubMed ID'], list):
                    experiment_object.pubmed_id = idf_dict['PubMed ID'][0]
                else:
                    experiment_object.pubmed_id = idf_dict['PubMed ID']
                experiment_object.has_publication = True

            experiment_object.save()

        platform_dict = {k: experiment[k]
                         for k in ('platform_accession_code', 'platform_accession_name')}
        return experiment_object, platform_dict

    def determine_sample_accession(self, experiment_accession: str, sample_source_name: str, sample_assay_name: str, filename: str) -> str:
        """Determine what to use as the sample's accession code.

        This is a complicated heuristic to determine the sample
        accession because there isn't a field that consistently
        contains it so we're trying to figure out a heuristic that
        will work for all the data. This may need even further
        refinements if we come across examples that break it.

        However, what's going on is that we think either the `source`
        or `assay` field will be the sample accession but it's not
        always the same.
        Ex: E-MEXP-669 has it in sample_assay_name.
        Therefore we try a few different things to determine which it
        is.

        The experiment accession must be prefixed since accessions
        are non-unique on AE, ex "Sample 1" is a valid assay name.

        """

        # It SEEMS like the filename often contains part or all of the
        # sample name so we first try to see if either field contains
        # the filename with the extension stripped off:
        if isinstance(filename, str):
            stripped_filename = ".".join(filename.split(".")[:-1])
            if stripped_filename != "":
                if stripped_filename in sample_source_name:
                    return experiment_accession + "-" + sample_source_name
                elif stripped_filename in sample_assay_name:
                    return experiment_accession + "-" + sample_assay_name

        # Accessions don't have spaces in them, but sometimes these
        # fields do so next we try to see if one has spaces and the
        # other doesn't:
        source_has_spaces = " " in sample_source_name
        assay_has_spaces = " " in sample_assay_name
        if assay_has_spaces and not source_has_spaces:
            return experiment_accession + "-" + sample_source_name
        elif source_has_spaces and not assay_has_spaces:
            return experiment_accession + "-" + sample_assay_name

        # We're out of options so return the longest one.
        if len(sample_source_name) >= len(sample_assay_name):
            return experiment_accession + "-" + sample_source_name
        else:
            return experiment_accession + "-" + sample_assay_name

    def create_samples_from_api(self,
                                experiment: Experiment,
                                platform_dict: Dict
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

        samples_endpoint = SAMPLES_URL.format(experiment.accession_code)
        r = utils.requests_retry_session().get(samples_endpoint, timeout=60)
        samples = r.json()["experiment"]["sample"]

        # The SDRF is the complete metadata record on a sample/property basis.
        # We run this through our harmonizer and then attach the properties
        # to our created samples.
        SDRF_URL_TEMPLATE = "https://www.ebi.ac.uk/arrayexpress/files/{code}/{code}.sdrf.txt"
        sdrf_url = SDRF_URL_TEMPLATE.format(code=experiment.accession_code)
        sdrf_samples = harmony.parse_sdrf(sdrf_url)
        harmonized_samples = harmony.harmonize(sdrf_samples)

        # An experiment can have many samples
        for sample in samples:

            # For some reason, this sample has no files associated with it.
            if "file" not in sample or len(sample['file']) == 0:
                continue

            # Each sample is given an experimenatlly-unique title.
            flat_sample = utils.flatten(sample)
            title = harmony.extract_title(flat_sample)

            # Figure out the Organism for this sample
            organism_name = UNKNOWN
            for characteristic in sample["characteristic"]:
                if characteristic["category"].upper() == "ORGANISM":
                    organism_name = characteristic["value"].upper()

            if organism_name == UNKNOWN:
                logger.warning("Sample from experiment %s did not specify the organism name.",
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

                # For ex: E-GEOD-15645
                if isinstance(sub_file['comment'], list):
                    sub_file_mod = sub_file
                    sub_file_mod['comment'] = sub_file['comment'][0]
                else:
                    sub_file_mod = sub_file

                # Some have the 'data' field, but not the actual data
                # Ex: E-GEOD-9656
                if sub_file_mod['type'] == "data" and sub_file_mod['comment'].get('value', None) != None:
                    has_raw = True
                if 'raw' in sub_file_mod['comment'].get('value', ''):
                    has_raw = True

            skip_sample = False
            for sub_file in sample['file']:

                # Skip derived data if we have it raw.
                if has_raw and "derived data" in sub_file['type']:
                    continue

                # XXX: This is a hack.
                # Don't get the raw data if it's only a 1-color sample.
                if 'Cy3' in str(sample) and 'Cy5' not in str(sample):
                    has_raw = False

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

                if not download_url:
                    logger.error("Sample %s did not specify a download url, skipping.",
                                 sample_accession_code,
                                 experiment_accession_code=experiment.accession_code,
                                 survey_job=self.survey_job.id)
                    skip_sample = True
                    continue

                if not filename:
                    logger.error("Sample %s did not specify a filename, skipping.",
                                 sample_accession_code,
                                 experiment_accession_code=experiment.accession_code,
                                 survey_job=self.survey_job.id)
                    skip_sample = True
                    continue

            if skip_sample:
                continue

            # The accession code is not a simple matter to determine.
            sample_source_name = sample["source"].get("name", "")
            sample_assay_name = sample["assay"].get("name", "")
            sample_accession_code = self.determine_sample_accession(
                experiment.accession_code,
                sample_source_name,
                sample_assay_name,
                filename)

            # Figure out the Organism for this sample
            organism_name = UNKNOWN
            for characteristic in sample["characteristic"]:
                if characteristic["category"].upper() == "ORGANISM":
                    organism_name = characteristic["value"].upper()

            if organism_name == UNKNOWN:
                logger.warning("Sample %s did not specify the organism name.",
                               sample_accession_code,
                               experiment_accession_code=experiment.accession_code,
                               survey_job=self.survey_job.id)
                organism = None
            else:
                organism = Organism.get_object_for_name(organism_name)

            # Create the sample object
            try:
                sample_object = Sample.objects.get(accession_code=sample_accession_code)
                logger.debug("Sample %s already exists, skipping object creation.",
                             sample_accession_code,
                             experiment_accession_code=experiment.accession_code,
                             survey_job=self.survey_job.id)
                continue
            except Sample.DoesNotExist:
                sample_object = Sample()

                # The basics
                sample_object.title = title
                sample_object.accession_code = sample_accession_code
                sample_object.source_archive_url = samples_endpoint
                sample_object.organism = organism
                sample_object.platform_name = platform_dict["platform_accession_name"]
                sample_object.platform_accession_code = platform_dict["platform_accession_code"]
                sample_object.technology = "MICROARRAY"
                sample_object.save()

                # Directly assign the harmonized properties
                harmonized_sample = harmonized_samples[title]
                for key, value in harmonized_sample.items():
                    setattr(sample_object, key, value)
                sample_object.save()

                sample_annotation = SampleAnnotation()
                sample_annotation.data = sample
                sample_annotation.sample = sample_object
                sample_annotation.is_ccdl = False
                sample_annotation.save()

                original_file = OriginalFile()
                original_file.source_filename = filename
                original_file.source_url = download_url
                original_file.is_downloaded = False
                original_file.is_archive = True
                original_file.has_raw = has_raw
                original_file.save()

                original_file_sample_association = OriginalFileSampleAssociation()
                original_file_sample_association.original_file = original_file
                original_file_sample_association.sample = sample_object
                original_file_sample_association.save()

            association = ExperimentSampleAssociation()
            association.experiment = experiment
            association.sample = sample_object
            association.save()

            # Create associations if they don't already exist
            try:
                assocation = ExperimentSampleAssociation.objects.get(
                    experiment=experiment, sample=sample_object)
            except ExperimentSampleAssociation.DoesNotExist:
                association = ExperimentSampleAssociation()
                association.experiment = experiment
                association.sample = sample_object
                association.save()

            try:
                assocation = ExperimentOrganismAssociation.objects.get(
                    experiment=experiment, organism=organism)
            except ExperimentOrganismAssociation.DoesNotExist:
                association = ExperimentOrganismAssociation()
                association.experiment = experiment
                association.organism = organism
                association.save()

            logger.info("Created Sample: " + str(sample_object),
                        experiment_accession_code=experiment.accession_code,
                        survey_job=self.survey_job.id)

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
            experiment, platform_dict = self.create_experiment_from_api(experiment_accession_code)
        except UnsupportedPlatformException as e:
            logger.info("Experiment was not on a supported platform, skipping.",
                        experiment_accession_code=experiment_accession_code,
                        survey_job=self.survey_job.id)
            return None, []

        samples = self.create_samples_from_api(experiment, platform_dict)
        return experiment, samples
