from datetime import datetime
from itertools import chain
from typing import Dict, List, Tuple

from data_refinery_common.enums import Downloaders
from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import (
    Experiment,
    ExperimentAnnotation,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    Sample,
    SampleAnnotation,
    SurveyJobKeyValue,
)
from data_refinery_common.utils import (
    get_normalized_platform,
    get_readable_affymetrix_names,
    get_supported_microarray_platforms,
)
from data_refinery_foreman.surveyor import harmony, utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor

logger = get_and_configure_logger(__name__)


EXPERIMENTS_URL = "https://www.ebi.ac.uk/biostudies/api/v1/studies/"
IDF_URL_TEMPLATE = "https://www.ebi.ac.uk/arrayexpress/files/{code}/{code}.idf.txt"
SDRF_URL_TEMPLATE = "https://www.ebi.ac.uk/biostudies/files/{code}/{code}.sdrf.txt"
SAMPLES_URL = EXPERIMENTS_URL + "{}/samples"
UNKNOWN = "UNKNOWN"


class UnsupportedPlatformException(Exception):
    pass


class ArrayExpressSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return Downloaders.ARRAY_EXPRESS.value

    def _parse_api_info_response(data: Dict) -> Dict:
        return {
            "file_count": data["files"],
            "ftp_link": data["ftpLink"],
            "last_updated_at": datetime.fromtimestamp(data["modified"] / 1000),
            "released_at": datetime.fromtimestamp(data["released"] / 1000),
        }

    @classmethod
    def _parse_api_data_response(cls, data: Dict) -> Dict:
        parsed_data = {}

        # Filter for supported section attributes.
        attribute_filter = lambda attr: attr["name"] in {
            "Description",
            "Organism",
            "Title",
        }
        attributes = data["section"]["attributes"]
        for attribute in filter(attribute_filter, attributes):
            parsed_data[attribute["name"].lower()] = attribute["value"]

        # Protocols.
        parsed_data["protocols"] = []
        protocol_filter = lambda entry: isinstance(entry, dict) and entry["type"] == "Protocols"
        value_defined_filter = lambda attr: "value" in attr
        subsections = data["section"]["subsections"]
        for subsection in filter(protocol_filter, chain.from_iterable(subsections)):
            # Protocols.
            protocol = {}
            for attribute in filter(value_defined_filter, subsection["attributes"]):
                protocol[attribute["name"].lower()] = attribute["value"]
            parsed_data["protocols"].append(protocol)

        # Array Designs.
        parsed_data["array_designs"] = []
        assay_filter = lambda entry: isinstance(entry, list) and entry["type"] == "Assays and Data"
        is_design = lambda entry: entry["type"] == "Array Designs"
        for subsection in filter(assay_filter, chain.from_iterable(subsections)):
            for assays_data_subsection in filter(is_design, subsection):
                links = chain.from_iterable(assays_data_subsection["links"])
                parsed_data["array_designs"].extend((link["url"] for link in links))

        return parsed_data

    @classmethod
    def _get_experiment_data(cls, accession_code: str) -> Dict:
        experiment_url = f"{EXPERIMENTS_URL}{accession_code}"
        data = {"url": experiment_url}

        # Get experiment main data.
        response = utils.requests_retry_session().get(experiment_url, timeout=60).json()
        data.update(cls._parse_api_data_response(response))

        # Get experiment info.
        response = (
            utils.requests_retry_session()
            .get(f"{EXPERIMENTS_URL}{accession_code}/info", timeout=60)
            .json()
        )
        info = cls._parse_api_info_response(response)

        data.update(info)
        return data

    @classmethod
    def _apply_metadata_to_experiment(cls, experiment: Experiment, metadata: Dict):
        experiment.title = metadata["title"]
        experiment.description = metadata.get("description", "Description not available.")
        experiment.source_database = "ARRAY_EXPRESS"
        # This will need to be updated if we ever use Array Express to get other kinds of data.
        experiment.technology = "MICROARRAY"
        experiment.source_first_published = metadata["released_at"]
        experiment.source_last_modified = metadata["last_updated_at"]

    def create_experiment_from_api(self, experiment_accession_code: str) -> Tuple[Experiment, Dict]:
        """Given an experiment accession code, create an Experiment object.

        Also returns a dictionary of additional information about the
        platform discovered for the experiment.

        Will raise an UnsupportedPlatformException if this experiment was
        conducted using a platform which we don't support.
        """

        try:
            experiment_data = self._get_experiment_data(accession_code=experiment_accession_code)
        except KeyError:
            logger.error(
                "Could not collect data from remote experiment source!",
                experiment_accession_code=experiment_accession_code,
                survey_job=self.survey_job.id,
            )
            raise

        array_designs = experiment_data["array_designs"]
        # This experiment has no platform at all, and is therefore useless.
        if not array_designs:
            logger.warn(
                "Remote experiment has no array design listed.",
                experiment_accession_code=experiment_accession_code,
                survey_job=self.survey_job.id,
            )
            raise UnsupportedPlatformException

        platforms = {}
        # If there is more than one array design listed in the experiment
        # then there is no other way to determine which array was used
        # for which sample other than looking at the header of the CEL
        # file. That obviously cannot happen until the CEL file has been
        # downloaded so we can just mark it as UNKNOWN and let the
        # downloader inspect the downloaded file to determine the
        # array then.
        if len(array_designs) != 1:
            platforms["platform_accession_code"] = UNKNOWN
            platforms["platform_accession_name"] = UNKNOWN
            platforms["manufacturer"] = UNKNOWN
        else:
            external_accession = array_designs[0]
            for platform in get_supported_microarray_platforms():
                if platform["external_accession"] == external_accession:
                    platforms["platform_accession_code"] = get_normalized_platform(
                        platform["platform_accession"]
                    )
                    # Illumina appears in the accession codes for
                    # platforms manufactured by Illumina
                    if "ILLUMINA" in platforms["platform_accession_code"].upper():
                        platforms["manufacturer"] = "ILLUMINA"
                        platforms["platform_accession_name"] = platform["platform_accession"]
                    # It's not Illumina, the only other supported Microarray platform is
                    # Affy. As our list of supported platforms grows this logic will
                    # need to get more sophisticated.
                    else:
                        platforms["manufacturer"] = "AFFYMETRIX"
                        platform_mapping = get_readable_affymetrix_names()
                        platforms["platform_accession_name"] = platform_mapping[
                            platform["platform_accession"]
                        ]

            if "platform_accession_code" not in platforms:
                # We don't know what platform this accession corresponds to.
                platforms["platform_accession_code"] = external_accession
                platforms["platform_accession_name"] = UNKNOWN
                platforms["manufacturer"] = UNKNOWN

        # Create the experiment.
        try:
            experiment = Experiment.objects.get(accession_code=experiment_accession_code)
            logger.debug(
                "Experiment already exists, skipping object creation.",
                experiment_accession_code=experiment_accession_code,
                survey_job=self.survey_job.id,
            )
        except Experiment.DoesNotExist:
            experiment = Experiment()
            experiment.accession_code = experiment_accession_code
            experiment.source_url = experiment_data["url"]
            self._apply_metadata_to_experiment(experiment, experiment_data)
            experiment.save()

            json_xa = ExperimentAnnotation()
            json_xa.experiment = experiment
            json_xa.data = experiment_data
            json_xa.is_ccdl = False
            json_xa.save()

            # Fetch and parse the IDF/SDRF file for any other fields.
            idf_url = IDF_URL_TEMPLATE.format(code=experiment_accession_code)
            idf_text = utils.requests_retry_session().get(idf_url, timeout=60).text

            idf_data = {}
            for line in idf_text.split("\n"):
                key_val = line.strip().split("\t")
                if len(key_val) == 2:
                    idf_data[key_val[0]] = key_val[1]
                elif len(key_val) > 2:
                    idf_data[key_val[0]] = key_val[1:]

            idf_xa = ExperimentAnnotation()
            idf_xa.data = idf_data
            idf_xa.experiment = experiment
            idf_xa.is_ccdl = False
            idf_xa.save()

            if "Investigation Title" in idf_data and isinstance(
                idf_data["Investigation Title"], str
            ):
                experiment.title = idf_data["Investigation Title"]

            if "Person Affiliation" in idf_data:
                # This is very rare, ex: E-MEXP-32
                if isinstance(idf_data["Person Affiliation"], list):
                    unique_people = list(set(idf_data["Person Affiliation"]))
                    experiment.submitter_institution = ", ".join(unique_people)[:255]
                else:
                    experiment.submitter_institution = idf_data["Person Affiliation"]

            if "Publication Title" in idf_data:
                # This will happen for some superseries.
                # Ex: E-GEOD-29536
                # Assume most recent is "best:, store the rest in experiment annotation.
                if isinstance(idf_data["Publication Title"], list):
                    experiment.publication_title = "; ".join(idf_data["Publication Title"])
                else:
                    experiment.publication_title = idf_data["Publication Title"]
                experiment.has_publication = True

            if "Publication DOI" in idf_data:
                experiment.has_publication = True
                if isinstance(idf_data["Publication DOI"], list):
                    experiment.publication_doi = ", ".join(idf_data["Publication DOI"])
                else:
                    experiment.publication_doi = idf_data["Publication DOI"]

            if "PubMed ID" in idf_data:
                experiment.has_publication = True
                if isinstance(idf_data["PubMed ID"], list):
                    experiment.pubmed_id = ", ".join(idf_data["PubMed ID"])
                else:
                    experiment.pubmed_id = idf_data["PubMed ID"]

            # Scrape publication title and authorship from Pubmed
            if experiment.pubmed_id:
                pubmed_metadata = utils.get_title_and_authors_for_pubmed_id(experiment.pubmed_id)
                experiment.publication_title = pubmed_metadata[0]
                experiment.publication_authors = pubmed_metadata[1]

            experiment.protocol_description = experiment_data["protocols"]
            experiment.save()

        return experiment, platforms

    def determine_sample_accession(
        self,
        experiment_accession: str,
        sample_source_name: str,
        sample_assay_name: str,
        filename: str,
    ) -> str:
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

    @staticmethod
    def extract_protocol_text(protocol_text):
        """Returns a string representation of protocol_text.

        protocol_text may be a string or a list containing both
        strings and dicts, like so (it's what the API returns
        sometimes, see E-MEXP-2381 as an example):
        [
          "Microarrays were imaged using an Agilent microarray scanner in XDR (eXtended Dynamic Range function) mode and a scan resolution of 5 \u00b5m.",
          {
            "br": null
          },
          "(Parameters: Scanning hardware = DNA Microarray Scanner BA [Agilent Technologies], Scanning software = Feature Extraction Software [Agilent])"
        ]
        """
        if not protocol_text:
            return ""
        elif type(protocol_text) == str:
            return protocol_text.strip()
        elif type(protocol_text) == list:
            # These can be {"br": None}, so skip non string lines
            return " ".join([line.strip() for line in protocol_text if type(line) == str])
        else:
            # Not sure what would get us here, but it's not worth raising an error over
            return str(protocol_text)

    @staticmethod
    def update_sample_protocol_info(existing_protocols, experiment_protocol, protocol_url):
        """Compares experiment_protocol with a sample's
        existing_protocols and updates the latter if the former includes
        any new entry.

        Returns a two-element tuple, the first is existing_protocols
        (which may or may not have been updated) and the second is a
        bool indicating whether exisiting_protocols has been updated.

        Note that the ArrayExpress experiment-level protocol may include
        multiple protocol entries.
        """

        if not "protocol" in experiment_protocol:
            return (existing_protocols, False)

        is_updated = False
        # Compare each entry in experiment protocol with the existing
        # protocols; if the entry is new, add it to exising_protocols.
        for new_protocol in experiment_protocol["protocol"]:
            new_protocol_text = new_protocol.get("text", "")
            new_protocol_text = ArrayExpressSurveyor.extract_protocol_text(new_protocol_text)

            # Ignore experiment-level protocols whose accession or text
            # field is unavailable or empty.
            if not new_protocol.get("accession", "").strip() or not new_protocol_text:
                continue

            new_protocol_is_found = False
            for existing_protocol in existing_protocols:
                if (
                    new_protocol.get("accession", "") == existing_protocol["Accession"]
                    and new_protocol_text == existing_protocol["Text"]
                    and new_protocol.get("type", "") == existing_protocol["Type"]
                ):
                    new_protocol_is_found = True
                    break
            if not new_protocol_is_found:
                existing_protocols.append(
                    {
                        "Accession": new_protocol["accession"],
                        "Text": new_protocol_text,
                        "Type": new_protocol.get("type", ""),  # in case 'type' field is unavailable
                        "Reference": protocol_url,
                    }
                )
                is_updated = True

        return (existing_protocols, is_updated)

    @staticmethod
    def _apply_harmonized_metadata_to_sample(sample: Sample, harmonized_metadata: dict):
        """Applies the harmonized metadata to `sample`"""
        for key, value in harmonized_metadata.items():
            setattr(sample, key, value)

    def create_samples_from_api(self, experiment: Experiment, platform_dict: Dict) -> List[Sample]:
        """Generates a Sample item for each sample in an AE experiment.

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

        title_field = harmony.determine_title_field(sdrf_samples, samples)
        harmonized_samples = harmony.harmonize_all_samples(sdrf_samples, title_field)

        # An experiment can have many samples
        for sample_data in samples:

            # For some reason, this sample has no files associated with it.
            if "file" not in sample_data or len(sample_data["file"]) == 0:
                continue

            # Each sample is given an experimenatlly-unique title.
            flat_sample = utils.flatten(sample_data)
            title = harmony.extract_title(flat_sample, title_field)

            # A sample may actually have many sub files.
            # If there is raw data, take that.
            # If not, take the derived.
            has_raw = False
            for sub_file in sample_data["file"]:

                # For ex: E-GEOD-15645
                if isinstance(sub_file["comment"], list):
                    sub_file_mod = sub_file
                    sub_file_mod["comment"] = sub_file["comment"][0]
                else:
                    sub_file_mod = sub_file

                # Some have the 'data' field, but not the actual data
                # Ex: E-GEOD-9656
                if (
                    sub_file_mod["type"] == "data"
                    and sub_file_mod["comment"].get("value", None) != None
                ):
                    has_raw = True

                # 'value' can be None, convert to an empty string to
                # make it easier to use.
                comment_value = sub_file_mod["comment"].get("value", "") or ""
                if "raw" in comment_value:
                    has_raw = True

            skip_sample = False
            for sub_file in sample_data["file"]:

                # Don't get the raw data if it's only a 1-color sample.
                if "Cy3" in str(sample_data) and "Cy5" not in str(sample_data):
                    has_raw = False

                # Skip derived data if we have it raw.
                if has_raw and "derived data" in sub_file["type"]:
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

                if not download_url:
                    logger.error(
                        "Sample %s did not specify a download url, skipping.",
                        sample_accession_code,
                        experiment_accession_code=experiment.accession_code,
                        survey_job=self.survey_job.id,
                        sub_file=sub_file,
                    )
                    skip_sample = True
                    continue

                if not filename:
                    logger.error(
                        "Sample %s did not specify a filename, skipping.",
                        sample_accession_code,
                        experiment_accession_code=experiment.accession_code,
                        survey_job=self.survey_job.id,
                        sub_file=sub_file,
                    )
                    skip_sample = True
                    continue

            if skip_sample:
                continue

            # The accession code is not a simple matter to determine.
            sample_source_name = sample_data["source"].get("name", "")
            sample_assay_name = sample_data["assay"].get("name", "")
            sample_accession_code = self.determine_sample_accession(
                experiment.accession_code,
                sample_source_name,
                sample_assay_name,
                filename,
            )

            # Figure out the Organism for this sample
            organism_name = UNKNOWN
            for characteristic in sample_data["characteristic"]:
                if characteristic["category"].upper() == "ORGANISM":
                    organism_name = characteristic["value"].upper()

            if organism_name == UNKNOWN:
                logger.error(
                    "Sample %s did not specify the organism name.",
                    sample_accession_code,
                    experiment_accession_code=experiment.accession_code,
                    survey_job=self.survey_job.id,
                )
                organism = None
                continue
            else:
                organism = Organism.get_object_for_name(organism_name)

            # Create the sample object
            try:
                # Associate it with the experiment, but since it
                # already exists it already has original files
                # associated with it and it's already been downloaded,
                # so don't add it to created_samples.
                sample_object = Sample.objects.get(accession_code=sample_accession_code)

                # If input experiment includes new protocol information,
                # update sample's protocol_info.
                existing_protocols = sample_object.protocol_info
                protocol_info, is_updated = self.update_sample_protocol_info(
                    existing_protocols,
                    experiment.protocol_description,
                    experiment.source_url + "/protocols",
                )
                if is_updated:
                    sample_object.protocol_info = protocol_info
                    sample_object.save()

                logger.debug(
                    "Sample %s already exists, skipping object creation.",
                    sample_accession_code,
                    experiment_accession_code=experiment.accession_code,
                    survey_job=self.survey_job.id,
                )
            except Sample.DoesNotExist:
                sample_object = Sample()

                # The basics
                sample_object.source_database = "ARRAY_EXPRESS"
                sample_object.title = title
                sample_object.accession_code = sample_accession_code
                sample_object.source_archive_url = samples_endpoint
                sample_object.organism = organism
                sample_object.platform_name = platform_dict["platform_accession_name"]
                sample_object.platform_accession_code = platform_dict["platform_accession_code"]
                sample_object.manufacturer = platform_dict["manufacturer"]
                sample_object.technology = "MICROARRAY"

                protocol_info, is_updated = self.update_sample_protocol_info(
                    existing_protocols=[],
                    experiment_protocol=experiment.protocol_description,
                    protocol_url=experiment.source_url + "/protocols",
                )
                # Do not check is_updated the first time because we must
                # save a list so we can append to it later.
                sample_object.protocol_info = protocol_info

                sample_object.save()

                # Directly assign the harmonized properties
                harmonized_sample = harmonized_samples[title]
                ArrayExpressSurveyor._apply_harmonized_metadata_to_sample(
                    sample_object, harmonized_sample
                )

                sample_annotation = SampleAnnotation()
                sample_annotation.name = "raw_metadata"
                sample_annotation.data = sample_data
                sample_annotation.sample = sample_object
                sample_annotation.is_ccdl = False
                sample_annotation.save()

                original_file = OriginalFile()
                original_file.filename = filename
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

                created_samples.append(sample_object)

                logger.debug(
                    "Created " + str(sample_object),
                    experiment_accession_code=experiment.accession_code,
                    survey_job=self.survey_job.id,
                    sample=sample_object.id,
                )

            # Create associations if they don't already exist
            ExperimentSampleAssociation.objects.get_or_create(
                experiment=experiment, sample=sample_object
            )

            ExperimentOrganismAssociation.objects.get_or_create(
                experiment=experiment, organism=organism
            )

        return created_samples

    def discover_experiment_and_samples(self) -> (Experiment, List[Sample]):

        experiment_accession_code = SurveyJobKeyValue.objects.get(
            survey_job_id=self.survey_job.id, key__exact="experiment_accession_code"
        ).value

        logger.info(
            "Surveying experiment with accession code: %s.",
            experiment_accession_code,
            survey_job=self.survey_job.id,
        )

        try:
            experiment, platform_dict = self.create_experiment_from_api(experiment_accession_code)
        except UnsupportedPlatformException:
            logger.info(
                "Experiment was not on a supported platform, skipping.",
                experiment_accession_code=experiment_accession_code,
                survey_job=self.survey_job.id,
            )
            return None, []
        except Exception:
            logger.exception(
                "Error occurred while surveying experiment!",
                experiment_accession_code=experiment_accession_code,
            )
            return None, []

        samples = self.create_samples_from_api(experiment, platform_dict)
        return experiment, samples
