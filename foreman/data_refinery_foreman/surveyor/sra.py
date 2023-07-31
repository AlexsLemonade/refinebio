import re
from typing import List
from urllib.parse import urlencode

from django.utils.dateparse import parse_date

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
    SurveyJob,
)
from data_refinery_common.rna_seq import _build_ena_file_url
from data_refinery_foreman.surveyor import harmony, utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor

logger = get_and_configure_logger(__name__)


DOWNLOAD_SOURCE = "NCBI"  # or "ENA". Change this to download from NCBI (US) or ENA (UK).
NCBI_DOWNLOAD_URL_TEMPLATE = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{accession}/{accession}"
ENA_URL_TEMPLATE = "https://www.ebi.ac.uk/ena/browser/view/{}"
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/"
ENA_BIOSAMPLE_API = "https://www.ebi.ac.uk/biosamples/samples"
ENA_PUBMED_API_TEMPLATE = "https://www.ebi.ac.uk/ena/xref/rest/json/search?accession={accession}&expanded=true&source=PubMed"


class UnsupportedDataTypeError(Exception):
    pass


class SraSurveyor(ExternalSourceSurveyor):

    """Surveys SRA for data.

    Implements the ExternalSourceSurveyor interface. It is worth
    noting that a large part of this class is parsing XML metadata
    about Experiments. The strategy for parsing the metadata was to take
    nearly all of the fields that could be extracted from the
    XML. Essentially all the XML parsing code is just working through
    the XML in the way it is formatted. For reference see the sample
    XML which has been included in ./test_sra_xml.py. This sample XML
    is used for testing, so the values that are extracted from it can
    be found in ./test_sra.py.
    """

    def source_type(self):
        return Downloaders.SRA.value

    @staticmethod
    def get_ena_json_api(endpoint: str, accession: str, params: dict = {}):

        params.update({"format": "json"})

        if "accession" not in params:
            params.update({"includeAccessions": accession})

        return f"{ENA_API}{endpoint}/?{urlencode(params, )}"

    @staticmethod
    def get_ena_json_api_response(endpoint: str, accession: str, params: dict = {}):
        response = utils.requests_retry_session().get(
            SraSurveyor.get_ena_json_api(endpoint, accession, params)
        )
        # handle logging of errors here consistently
        if response.status_code != 200:
            logger.error(
                f"Could not discover {endpoint}.",
                accession=accession,
                status_code=response.status_code,
                url=response.url,
            )
            return {}
        return response.json()

    @staticmethod
    def gather_experiment_metadata(experiment_accession: str):
        experiments = SraSurveyor.get_ena_json_api_response(
            "search",
            experiment_accession,
            {
                "result": "study",
                "fields": "all",
            },
        )
        if len(experiments) > 0:
            return experiments[0]
        return {}

    @staticmethod
    def gather_experiment_runs_metadata(accession: str):
        # if accession is experiment it will return all samples
        # if accession is sample it will return just the sample
        samples_metadata = SraSurveyor.get_ena_json_api_response(
            "search",
            accession,
            {
                "result": "read_run",
                "fields": "all",
            },
        )

        samples_metadata = [SraSurveyor.gather_biosample_metadata(s) for s in samples_metadata]

        return samples_metadata

    @staticmethod
    def gather_file_report(run_accession: str):
        return SraSurveyor.get_ena_json_api_response(
            "filereport",
            run_accession,
            {
                "accession": run_accession,
                "result": "read_run",
                "fields": "all",
            },
        )

    @staticmethod
    def gather_experiment_biosample_metadata(biosample_accession: str):
        response = utils.requests_retry_session().get(
            f"{ENA_BIOSAMPLE_API}/{biosample_accession}?format=json"
        )

        if response.status_code != 200:
            logger.error(
                "Could not discover Biosample.",
                accession=biosample_accession,
                status_code=response.status_code,
                url=response.url,
            )
            return {}
        return response.json()

    @staticmethod
    def gather_experiment_pubmed_metadata(accession: str):
        response = utils.requests_retry_session().get(
            ENA_PUBMED_API_TEMPLATE.format(accession=accession)
        )

        if response.status_code != 200:
            logger.error(
                "Could not discover Pubmed ID.",
                accession=accession,
                status_code=response.status_code,
                url=response.url,
            )
            return []
        return response.json()

    @staticmethod
    def gather_biosample_metadata(sample_metadata):
        biosample_metadata = SraSurveyor.gather_experiment_biosample_metadata(
            sample_metadata["sample_accession"]
        )
        biosample_metadata.pop(
            "accession", None
        )  # dont overwrite accession with biosample accession
        biosample_metadata.pop("status", None)  # these should match but biosample is all uppercase
        return {**sample_metadata, **utils.flatten_dict(biosample_metadata)}

    @staticmethod
    def gather_all_metadata(accession: str):
        # check if experiment, if not fetch the sample, biosample and experiment
        samples_metadata = SraSurveyor.gather_experiment_runs_metadata(accession)

        if not re.match(r"^(S|E|D)RP\d+", accession):
            experiment_accession = samples_metadata[0]["secondary_study_accession"]
        else:
            experiment_accession = accession

        experiment_metadata = SraSurveyor.gather_experiment_metadata(experiment_accession)

        if "accession" in experiment_metadata:
            pubmed_metadata = SraSurveyor.gather_experiment_pubmed_metadata(
                experiment_metadata["accession"]
            )
            if len(pubmed_metadata) > 0:
                experiment_metadata["pubmed_id"] = pubmed_metadata[0]["Source Primary Accession"]

        for sample_metadata in samples_metadata:
            sample_metadata.update(utils.flatten_dict(experiment_metadata, "study_ena"))

        return experiment_metadata, samples_metadata

    @staticmethod
    def _build_ncbi_file_url(run_accession: str):
        """Build the path to the hypothetical .sra file we want"""
        return NCBI_DOWNLOAD_URL_TEMPLATE.format(accession=run_accession)

    @staticmethod
    def _apply_harmonized_metadata_to_sample(sample: Sample, metadata: dict):
        """Harmonizes the metadata and applies it to `sample`"""
        harmonizer = harmony.Harmonizer()
        harmonized_sample = harmonizer.harmonize_sample(metadata)
        for key, value in harmonized_sample.items():
            setattr(sample, key, value)

    @staticmethod
    def _apply_metadata_to_experiment(experiment: Experiment, metadata: dict):
        experiment.source_url = ENA_URL_TEMPLATE.format(experiment.accession_code)
        experiment.source_database = "SRA"
        experiment.technology = "RNA-SEQ"

        experiment.title = metadata.get("study_title", "")
        experiment.description = utils.get_nonempty(
            metadata, "study_description", "No description."
        )
        experiment.submitter_institution = utils.get_nonempty(metadata, "center_name", "")
        experiment.source_first_published = parse_date(metadata["first_public"])
        experiment.source_last_modified = parse_date(metadata["last_updated"])
        experiment.alternate_accession_code = utils.get_nonempty(metadata, "geo_accession", None)
        experiment.protocol_description = utils.get_nonempty(metadata, "protocol", [])

        experiment.pubmed_id = metadata.get("pubmed_id", "")
        experiment.has_publication = "pubmed_id" in metadata
        # Scrape publication title and authorship from Pubmed
        if experiment.pubmed_id:
            title, authors = utils.get_title_and_authors_for_pubmed_id(experiment.pubmed_id)
            experiment.publication_title = title
            experiment.publication_authors = authors

    @staticmethod
    def get_files_urls(run_accession, is_paired):
        if DOWNLOAD_SOURCE == "NCBI":
            return [SraSurveyor._build_ncbi_file_url(run_accession)]
        if is_paired:
            return [
                _build_ena_file_url(run_accession, "_1"),
                _build_ena_file_url(run_accession, "_2"),
            ]
        return [_build_ena_file_url(run_accession)]

    @staticmethod
    def _generate_original_files(sample, sample_metadata):
        is_paired = sample_metadata["library_layout"] == "PAIRED"
        files_urls = SraSurveyor.get_files_urls(sample.accession_code, is_paired)
        file_reports = SraSurveyor.gather_file_report(sample.accession_code)

        if len(file_reports) <= 2:
            for file_url in files_urls:
                original_file, _ = OriginalFile.objects.get_or_create(
                    source_url=file_url, source_filename=file_url.split("/")[-1], has_raw=True
                )
                OriginalFileSampleAssociation.objects.get_or_create(
                    original_file=original_file, sample=sample
                )
        else:
            # If there's 3 fastq files in ENA's FTP server that's
            # because one is the unmated reads and we can't
            # reliably convert from the .SRA. Therefore set the
            # file_urls to be the ENA FTP URLs.
            for file_report in file_reports:
                file_url = file_report["fastq_ftp"]

                # Skip the one for unmated reads: the unmated
                # reads file lacks a _X postfix. It's be better if
                # it was clearly labeled, but it seems to be
                # consistent. The unmated reads file also seems to
                # be the smallest by far, but I'm unsure how
                # reliable that is.
                if "_1.fastq.gz" in file_url or "_2.fastq.gz" in file_url:
                    original_file, _ = OriginalFile.objects.get_or_create(
                        source_url=file_url,
                        source_filename=file_url.split("/")[-1],
                        has_raw=True,
                        expected_size_in_bytes=file_report["fastq_bytes"],
                        expected_md5=file_report["fastq_md5"],
                    )
                    OriginalFileSampleAssociation.objects.get_or_create(
                        original_file=original_file, sample=sample
                    )

    def _generate_experiment_and_samples(
        self, experiment_or_sample_accession: str
    ) -> (Experiment, List[Sample]):
        """Generates Experiments and Samples for the provided run_accession."""
        experiment_metadata, samples_metadata = SraSurveyor.gather_all_metadata(
            experiment_or_sample_accession
        )

        experiment = self._generate_experiment(experiment_metadata)
        samples = [
            self._generate_sample(experiment, sample_metadata)
            for sample_metadata in samples_metadata
        ]

        return experiment, samples

    def _generate_experiment(self, experiment_metadata: dict) -> Experiment:
        """Generates an Experiment for the provided experiment_accession."""

        if not experiment_metadata:
            return None

        # secondary_study_accession is DRP, SRP etc
        # study_accession is PRJDB etc
        experiment_accession = experiment_metadata.get("secondary_study_accession")

        experiment, is_created = Experiment.objects.get_or_create(
            accession_code=experiment_accession
        )

        if is_created:
            SraSurveyor._apply_metadata_to_experiment(experiment, experiment_metadata)
            experiment.save()
            experiment_annotation = ExperimentAnnotation()
            experiment_annotation.experiment = experiment
            experiment_annotation.data = experiment_metadata
            experiment_annotation.is_ccdl = False
            experiment_annotation.save()
        else:
            logger.debug(
                "Experiment already exists, skipping object creation.",
                experiment_accession_code=experiment_accession,
                survey_job=self.survey_job.id,
            )

        return experiment

    def _generate_sample(self, experiment: Experiment, sample_metadata: dict) -> Sample:
        """Generates or updates a Sample for the provided Experiment and sample metadata."""

        run_accession = sample_metadata["run_accession"]

        # Figure out the Organism for this sample
        organism_name = sample_metadata.pop("scientific_name", None)
        if not organism_name:
            logger.error("Could not discover organism type for run.", accession=run_accession)
            return (None, None)  # This will cascade properly

        organism = Organism.get_object_for_name(organism_name)

        sample_accession = sample_metadata.pop("run_accession")
        sample, sample_created = Sample.objects.get_or_create(accession_code=sample_accession)

        if sample_created:
            sample.source_database = "SRA"
            sample.organism = organism
            sample.manufacturer = sample_metadata.get("instrument_platform", "UNKNOWN")
            sample.platform_name = sample_metadata.get("instrument_model", "UNKNOWN")
            sample.platform_accession_code = sample.platform_name.replace(" ", "")
            sample.technology = "RNA-SEQ"

            # In the rare case that we get something we don't expect
            # we don't want to process this sample so we can set it to UNKNOWN.
            if sample.manufacturer not in ["ILLUMINA", "ION_TORRENT"]:
                sample.manufacturer = "UNKNOWN"
                sample.platform_name = "UNKNOWN"
                sample.platform_accession_code = "UNKNOWN"
                logger.debug(
                    "Sample %s has unrecognized manufacturer.",
                    sample_accession_code=sample_accession,
                    experiment_accession_code=experiment.accession_code,
                    survey_job=self.survey_job.id,
                )

            SraSurveyor._apply_harmonized_metadata_to_sample(sample, sample_metadata)

            self._generate_original_files(sample, sample_metadata)

            sample.save()
        else:
            logger.debug(
                "Sample %s already exists, skipping object creation.",
                sample_accession_code=sample_accession,
                experiment_accession_code=experiment.accession_code,
                survey_job=self.survey_job.id,
            )

        # i think we should be collecting the sample protocols
        # on the experiment here instead
        protocol_info, protocol_updated = self.update_sample_protocol_info(
            existing_protocols=sample.protocol_info,
            experiment_protocol=experiment.protocol_description,
            experiment_url=experiment.source_url,
        )

        if sample_created or protocol_updated:
            sample.protocol_info = protocol_info
            sample.save()

        # Create associations if they don't already exist
        ExperimentSampleAssociation.objects.get_or_create(experiment=experiment, sample=sample)

        ExperimentOrganismAssociation.objects.get_or_create(
            experiment=experiment, organism=organism
        )

        return sample

    @staticmethod
    def update_sample_protocol_info(existing_protocols, experiment_protocol, experiment_url):
        """Compares experiment_protocol with a sample's
        existing_protocols and update the latter if the former is new.

        Returns a tuple whose first element is existing_protocols (which
        may or may not have been updated) and the second is a boolean
        indicating whether exisiting_protocols has been updated.
        """
        # ensure we are working with a list
        # since the default value for a JSON field is an empty dict
        existing_protocols_list = []
        if existing_protocols is not {}:
            existing_protocols_list = existing_protocols

        if experiment_protocol == "Protocol was never provided.":
            return (existing_protocols_list, False)

        existing_descriptions = [protocol["Description"] for protocol in existing_protocols_list]
        if experiment_protocol in existing_descriptions:
            return (existing_protocols_list, False)

        existing_protocols_list.append(
            {"Description": experiment_protocol, "Reference": experiment_url}
        )
        return (existing_protocols_list, True)

    def discover_experiment_and_samples(self):
        """Returns an experiment and a list of samples for an SRA accession"""
        survey_job = SurveyJob.objects.get(id=self.survey_job.id)
        survey_job_properties = survey_job.get_properties()
        accession = survey_job_properties["experiment_accession_code"]

        logger.debug("Surveying SRA Accession %s", accession, survey_job=self.survey_job.id)
        return self._generate_experiment_and_samples(accession)
