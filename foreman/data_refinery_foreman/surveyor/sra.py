import requests
from typing import List, Dict
import xml.etree.ElementTree as ET
from data_refinery_common.models import (
    Batch,
    BatchKeyValue,
    File,
    SurveyJob,
    SurveyJobKeyValue,
    Organism
)
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger


logger = get_and_configure_logger(__name__)


DDBJ_URL_BASE = "ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/"
ENA_URL_TEMPLATE = "https://www.ebi.ac.uk/ena/data/view/{}&display=xml"
NCBI_DOWNLOAD_URL_TEMPLATE = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=dload&run_list={}&format=fastq"  # noqa
TEST_XML = "SRA200/SRA200001"
TAIL = "SRA200001.run.xml"

# SRA has a convention around its accession IDs. Each accession ID has
# a three letter prefix. The first is what organization the data was
# submitted to:
# A - America
# E - Europe
# D - Japan

# The second letter is presumably related to the type of data, but for
# RNA data it is always R so I haven't needed to verify this. The
# third letter corresponds to what type of data the accession ID is
# referring to:
SUBMISSION_CODE = "A"
EXPERIMENT_CODE = "X"
STUDY_CODE = "P"
SAMPLE_CODE = "S"
RUN_CODE = "R"


class UnsupportedDataTypeError(BaseException):
    pass


class SraSurveyor(ExternalSourceSurveyor):
    def source_type(self):
        return Downloaders.SRA.value

    def determine_pipeline(self,
                           batch: Batch,
                           key_values: List[BatchKeyValue] = []):
        return ProcessorPipeline.SALMON

    @staticmethod
    def gather_submission_metadata(metadata: Dict) -> None:
        response = requests.get(ENA_URL_TEMPLATE.format(metadata["submission_accession"]))
        submission_xml = ET.fromstring(response.text)[0]
        submission_metadata = submission_xml.attrib

        # We already have these
        submission_metadata.pop("accession")
        submission_metadata.pop("alias")

        metadata.update(submission_metadata)

        for child in submission_xml:
            if child.tag == "TITLE":
                metadata["submission_title"] = child.text
            elif child.tag == "SUBMISSION_ATTRIBUTES":
                for grandchild in child:
                    metadata[grandchild.find("TAG").text.lower()] = grandchild.find("VALUE").text

    @staticmethod
    def gather_library_metadata(metadata: Dict, library: ET.Element) -> None:
        for child in library:
            if child.tag == "LIBRARY_LAYOUT":
                metadata["library_layout"] = child[0].tag
            else:
                metadata[child.tag.lower()] = child.text

        # SRA contains several types of data. We only want RNA-Seq for now.
        if metadata["library_strategy"] != "RNA-Seq":
            raise UnsupportedDataTypeError("library_strategy not RNA-Seq.")
        if metadata["library_source"] not in ["TRANSCRIPTOMIC", "OTHER"]:
            raise UnsupportedDataTypeError("library_source not TRANSCRIPTOMIC or OTHER.")

    @staticmethod
    def parse_read_spec(metadata: Dict, read_spec: ET.Element, counter: int) -> None:
        for child in read_spec:
            key = "read_spec_{}_{}".format(str(counter), child.tag.replace("READ_", "").lower())
            metadata[key] = child.text

    @staticmethod
    def gather_spot_metadata(metadata: Dict, spot: ET.Element) -> None:
        """We don't actually know what this metadata is about."""
        for child in spot:
            if child.tag == "SPOT_DECODE_SPEC":
                read_spec_counter = 0
                for grandchild in child:
                    if grandchild.tag == "SPOT_LENGTH":
                        metadata["spot_length"] = grandchild.text
                    elif grandchild.tag == "READ_SPEC":
                        SraSurveyor.parse_read_spec(metadata, grandchild, read_spec_counter)
                        read_spec_counter = read_spec_counter + 1

    @staticmethod
    def gather_experiment_metadata(metadata: Dict) -> None:
        response = requests.get(ENA_URL_TEMPLATE.format(metadata["experiment_accession"]))
        experiment_xml = ET.fromstring(response.text)

        experiment = experiment_xml[0]
        for child in experiment:
            if child.tag == "TITLE":
                metadata["experiment_title"] = child.text
            elif child.tag == "DESIGN":
                for grandchild in child:
                    if grandchild.tag == "DESIGN_DESCRIPTION":
                        metadata["experiment_desing_description"] = grandchild.text
                    elif grandchild.tag == "LIBRARY_DESCRIPTOR":
                        SraSurveyor.gather_library_metadata(metadata, grandchild)
                    elif grandchild.tag == "SPOT_DESCRIPTOR":
                        SraSurveyor.gather_spot_metadata(metadata, grandchild)
            elif child.tag == "PLATFORM":
                # This structure is extraneously nested.
                # This is used as the platform_accession_code for SRA
                # batches, which becomes part of file paths, so we
                # don't want any spaces in it.
                metadata["platform_instrument_model"] = child[0][0].text.replace(" ", "")

    @staticmethod
    def parse_run_link(run_link: ET.ElementTree) -> (str, str):
        key = ""
        value = ""

        # The first level in this element is a XREF_LINK which is
        # really just an unnecessary level
        for child in run_link[0]:
            if child.tag == "DB":
                key = child.text.lower().replace("ena-", "") + "_accession"
            elif child.tag == "ID":
                value = child.text

        return (key, value)

    @staticmethod
    def parse_attribute(attribute: ET.ElementTree, key_prefix: str ="") -> (str, str):
        key = ""
        value = ""

        for child in attribute:
            if child.tag == "TAG":
                key = key_prefix + child.text.lower().replace("-", "_").replace(" ", "_")
            elif child.tag == "VALUE":
                value = child.text

        return (key, value)

    @staticmethod
    def gather_run_metadata(run_accession: str) -> Dict:
        """A run refers to a specific read in an experiment."""
        discoverable_accessions = ["study_accession", "sample_accession", "submission_accession"]
        response = requests.get(ENA_URL_TEMPLATE.format(run_accession))
        run_xml = ET.fromstring(response.text)
        run = run_xml[0]

        useful_attributes = ["center_name", "run_center", "run_date", "broker_name"]
        metadata = {key: run.attrib[key] for key in useful_attributes}
        metadata["run_accession"] = run_accession

        for child in run:
            if child.tag == "EXPERIMENT_REF":
                metadata["experiment_accession"] = child.attrib["accession"]
            elif child.tag == "RUN_LINKS":
                for grandchild in child:
                    key, value = SraSurveyor.parse_run_link(grandchild)
                    if value != "" and key in discoverable_accessions:
                        metadata[key] = value
            elif child.tag == "RUN_ATTRIBUTES":
                for grandchild in child:
                    key, value = SraSurveyor.parse_attribute(grandchild, "run_")
                    metadata[key] = value

        return metadata

    @staticmethod
    def gather_sample_metadata(metadata: Dict) -> None:
        response = requests.get(ENA_URL_TEMPLATE.format(metadata["sample_accession"]))
        sample_xml = ET.fromstring(response.text)

        sample = sample_xml[0]
        metadata["sample_center_name"] = sample.attrib["center_name"]
        for child in sample:
            if child.tag == "TITLE":
                metadata["sample_title"] = child.text
            elif child.tag == "SAMPLE_NAME":
                for grandchild in child:
                    if grandchild.tag == "TAXON_ID":
                        metadata["organism_id"] = grandchild.text
                    elif grandchild.tag == "SCIENTIFIC_NAME":
                        metadata["organism_name"] = grandchild.text.upper()
            elif child.tag == "SAMPLE_ATTRIBUTES":
                for grandchild in child:
                    key, value = SraSurveyor.parse_attribute(grandchild, "sample_")
                    metadata[key] = value

    @staticmethod
    def gather_study_metadata(metadata: Dict) -> None:
        metadata = {"study_accession": "DRP000595"}
        response = requests.get(ENA_URL_TEMPLATE.format(metadata["study_accession"]))
        study_xml = ET.fromstring(response.text)

        study = study_xml[0]
        for child in study:
            if child.tag == "DESCRIPTOR":
                for grandchild in child:
                    # STUDY_TYPE is the only tag which uses attributes
                    # instead of the text for whatever reason
                    if grandchild.tag == "STUDY_TYPE":
                        metadata[grandchild.tag.lower()] = grandchild.attrib["existing_study_type"]
                    else:
                        metadata[grandchild.tag.lower()] = grandchild.text
            elif child.tag == "STUDY_ATTRIBUTES":
                for grandchild in child:
                    key, value = SraSurveyor.parse_attribute(grandchild, "study_")
                    metadata[key] = value

    @staticmethod
    def gather_all_metadata(run_accession):
        metadata = SraSurveyor.gather_run_metadata(run_accession)

        SraSurveyor.gather_experiment_metadata(metadata)
        SraSurveyor.gather_sample_metadata(metadata)
        SraSurveyor.gather_study_metadata(metadata)
        SraSurveyor.gather_submission_metadata(metadata)

        return metadata

    def _generate_batch(self, run_accession: str) -> None:
        """Generates a Batch for each sample in samples.

        Uses the metadata contained in experiment (which should be
        generated via get_experiment_metadata) to add additional
        metadata to each Batch. If replicate_raw is True (the default)
        then only raw files will be replicated. Otherwise all files
        will be replicated.
        """
        metadata = SraSurveyor.gather_all_metadata(run_accession)

        file = File(name=(run_accession + ".fastq.gz"),
                    download_url=NCBI_DOWNLOAD_URL_TEMPLATE.format(run_accession),
                    raw_format=".fastq.gz",
                    processed_format=".tar.gz",
                    size_in_bytes=-1)  # Will have to be determined later

        self.add_batch(platform_accession_code=metadata.pop("platform_instrument_model"),
                       experiment_accession_code=metadata.pop("experiment_accession"),
                       organism_id=metadata.pop("organism_id"),
                       organism_name=metadata.pop("organism_name"),
                       experiment_title=metadata.pop("experiment_title"),
                       release_date=metadata.pop("run_ena_first_public"),
                       last_uploaded_date=metadata.pop("run_ena_last_update"),
                       files=[file],
                       key_values=metadata)

    @staticmethod
    def get_next_accession(last_accession: str) -> str:
        """Increments a SRA accession number by one.

        E.g. if last_accession is "DRR002116" then "DRR002117" will be
        returned.
        """
        prefix = last_accession[0:3]
        number = int(last_accession[3:])
        number = number + 1
        return prefix + "{0:06d}".format(number)

    def discover_batches(self):
        survey_job = SurveyJob.objects.get(id=self.survey_job.id)
        survey_job_properties = survey_job.get_properties()

        logger.info("Surveying SRA runs with accessions in the range of %s to %s.",
                    survey_job_properties["start_accession"],
                    survey_job_properties["end_accession"])

        current_accession = survey_job_properties["start_accession"]

        # By evaluating this conditional at the end of the loop
        # instead of the beginning, we achieve the functionality of a
        # do-while loop.
        surveyed_last_accession = False
        if not surveyed_last_accession:
            # Change to debug
            logger.info("Surveying SRA Run Accession %d", survey_job=survey_job.id)
            try:
                self._generate_batch(current_accession)
            except Exception as e:
                logger.exception("Exception caught while trying to generate a batch.",
                                 survey_job=self.survey_job.id,
                                 run_accession=current_accession)
                return False

            surveyed_last_accession = current_accession == survey_job_properties["end_accession"]
            current_accession = SraSurveyor.get_next_accession(current_accession)

        return True
