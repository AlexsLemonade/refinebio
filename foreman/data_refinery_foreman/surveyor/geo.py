import dateutil.parser
import GEOparse
import requests

from re import sub, split, match
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
    OriginalFileSampleAssociation
)
from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

class GeoUnsupportedPlatformException(BaseException):
    pass

class GeoSurveyor(ExternalSourceSurveyor):
    """Surveys NCBI GEO for data.

    Implements the GEO interface.
    """

    def source_type(self):
        return Downloaders.GEO.value

    def get_raw_url(self, experiment_accession_code):
        """ """
        geo = experiment_accession_code.upper()
        geotype = geo[:3]
        range_subdir = sub(r"\d{1,3}$", "nnn", geo)

        raw_url_template = ("ftp://ftp.ncbi.nlm.nih.gov/geo/"
                  "{root}/{range_subdir}/{record}/suppl/{record_file}")
        raw_url = raw_url_template.format(root="series",
                            range_subdir=range_subdir,
                            record=geo,
                            record_file="%s_RAW.tar" % geo)

        return raw_url

    def get_non_normalized_url(self, experiment_accession_code):
        """ """
        geo = experiment_accession_code.upper()
        geotype = geo[:3]
        range_subdir = sub(r"\d{1,3}$", "nnn", geo)

        nn_url_template = ("ftp://ftp.ncbi.nlm.nih.gov/geo/"
                  "{root}/{range_subdir}/{record}/suppl/{record_file}")
        nn_url = nn_url_template.format(root="series",
                            range_subdir=range_subdir,
                            record=geo,
                            record_file="%s_non-normalized.txt.gz" % geo)

        return nn_url

    def get_miniml_url(self, experiment_accession_code):
        """ """
        geo = experiment_accession_code.upper()
        geotype = geo[:3]
        range_subdir = sub(r"\d{1,3}$", "nnn", geo)

        min_url_template = ("ftp://ftp.ncbi.nlm.nih.gov/geo/"
                  "{root}/{range_subdir}/{record}/miniml/{record_file}")
        min_url = min_url_template.format(root="series",
                            range_subdir=range_subdir,
                            record=geo,
                            record_file="%s_family.xml.tgz" % geo)

        return min_url

    def create_experiment_and_samples_from_api(self, experiment_accession_code) -> (Experiment, List[Sample]):
        """ """

        # XXX: Maybe we should have an EFS tmp? This could potentially fill up if not tracked.
        # Cleaning up is tracked here: https://github.com/guma44/GEOparse/issues/41
        gse = GEOparse.get_GEO(experiment_accession_code, destdir='/tmp')

        # Create the experiment object
        try:
            experiment_object = Experiment.objects.get(accession_code=experiment_accession_code)
            logger.error("Experiment %s already exists, skipping object creation.",
                experiment_accession_code,
                survey_job=self.survey_job.id)
        except Experiment.DoesNotExist:
            experiment_object = Experiment()
            experiment_object.accession_code = experiment_accession_code
            experiment_object.source_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + experiment_accession_code
            experiment_object.source_database = "GEO"
            experiment_object.name = gse.metadata.get('title', [''])[0]
            experiment_object.description = gse.metadata.get('summary', [''])[0]
            experiment_object.platform_name = gse.metadata["platform_id"][0] # TODO: Lookup GEO-GPL
            experiment_object.platform_accession_code = gse.metadata["platform_id"][0]
            experiment_object.source_first_published = dateutil.parser.parse(gse.metadata["submission_date"][0] + " 00:00:00 UTC")
            experiment_object.source_last_updated = dateutil.parser.parse(gse.metadata["last_update_date"][0] + " 00:00:00 UTC")
            experiment_object.submitter_institution = ", ".join(list(set(gse.metadata["contact_institute"])))
            experiment_object.pubmed_id = gse.metadata.get("pubmed_id", [""])[0]
            experiment_object.save()

            experiment_annotation = ExperimentAnnotation()
            experiment_annotation.data = gse.metadata
            experiment_annotation.experiment = experiment_object
            experiment_annotation.is_ccdl = False
            experiment_annotation.save()

        has_raw = True

        # Okay, here's the situation!
        # Sometimes, samples have a direct single representation for themselves.
        # Othertimes, there is a single file with references to every sample in it.

        all_samples = []
        for sample_accession_code, sample in gse.gsms.items():

            try:
                sample_object = Sample.objects.get(accession_code=sample_accession_code)
                logger.error("Sample %s from experiment %s already exists, skipping object creation.",
                         sample_accession_code,
                         experiment_object.accession_code,
                         survey_job=self.survey_job.id)

                all_samples.append(sample_object)

                association = ExperimentSampleAssociation()
                association.experiment = experiment_object
                association.sample = sample_object
                association.save()
                continue
            except Sample.DoesNotExist:
                organism = Organism.get_object_for_name(sample.metadata['organism_ch1'][0].upper())

                # TODO: This is incomplete
                sample_object = Sample()
                sample_object.accession_code = sample_accession_code
                sample_object.organism = organism
                sample_object.title = sample.metadata['title'][0]
                sample_object.save()

                all_samples.append(sample_object)

                logger.info("Created Sample: " + str(sample_object))

                sample_annotation = SampleAnnotation()
                sample_annotation.sample = sample_object
                sample_annotation.data = sample.metadata
                sample_annotation.is_ccdl = False
                sample_annotation.save()

                sample_supplements = sample.metadata.get('supplementary_file', [])
                for supplementary_file_url in sample_supplements:

                    # Why do they give us this?
                    if supplementary_file_url == "NONE":
                        break

                    # We never want these!
                    if "idat.gz" in supplementary_file_url:
                        continue

                    try:
                        original_file = OriginalFile.objects.get(source_url=supplementary_file_url)
                    except OriginalFile.DoesNotExist:
                        original_file = OriginalFile()
                        # So - this is _usually_ true, but not always. I think it's submitter supplied.
                        original_file.source_filename = supplementary_file_url.split('/')[-1]
                        original_file.source_url = supplementary_file_url
                        original_file.is_downloaded = False
                        original_file.is_archive = True
                        original_file.has_raw = has_raw
                        original_file.save()

                        logger.info("Created OriginalFile: " + str(original_file))

                    original_file_sample_association = OriginalFileSampleAssociation()
                    original_file_sample_association.sample = sample_object
                    original_file_sample_association.original_file = original_file
                    original_file_sample_association.save()

                try:
                    assocation = ExperimentSampleAssociation.objects.get(experiment=experiment_object, sample=sample_object)
                except ExperimentSampleAssociation.DoesNotExist:
                    association = ExperimentSampleAssociation()
                    association.experiment = experiment_object
                    association.sample = sample_object
                    association.save()

        # These may or may not contain what we want.
        for experiment_supplement_url in gse.metadata.get('supplementary_file', []):

            try:
                original_file = OriginalFile.objects.get(source_url=experiment_supplement_url)
            except OriginalFile.DoesNotExist:
                original_file = OriginalFile()

                # So - this is _usually_ true, but not always. I think it's submitter supplied.
                original_file.source_filename = experiment_supplement_url.split('/')[-1]
                original_file.source_url = experiment_supplement_url
                original_file.is_downloaded = False
                original_file.is_archive = True
                original_file.has_raw = has_raw
                original_file.save()

                logger.info("Created OriginalFile: " + str(original_file))

            for sample_object in all_samples:
                original_file_sample_association = OriginalFileSampleAssociation()
                original_file_sample_association.sample = sample_object
                original_file_sample_association.original_file = original_file
                original_file_sample_association.save()

        # These are the Miniml/Soft/Matrix URLs that are always(?) provided?
        for family_url in [self.get_miniml_url(experiment_accession_code)]:

            try:
                original_file = OriginalFile.objects.get(source_url=family_url)
            except OriginalFile.DoesNotExist:
                original_file = OriginalFile()
                original_file.source_filename = family_url.split('/')[-1]
                original_file.source_url = family_url
                original_file.is_downloaded = False
                original_file.is_archive = True
                original_file.has_raw = has_raw
                original_file.save()
                logger.info("Created OriginalFile: " + str(original_file))

            for sample_object in all_samples:
                original_file_sample_association = OriginalFileSampleAssociation()
                original_file_sample_association.sample = sample_object
                original_file_sample_association.original_file = original_file
                original_file_sample_association.save()

        return experiment_object, all_samples

    def discover_experiment_and_samples(self) -> (Experiment, List[Sample]):
        """

        """

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
            experiment, samples = self.create_experiment_and_samples_from_api(experiment_accession_code)
        except GeoUnsupportedPlatformException as e:
            logger.info("Experiment with accession code: %s was not on a supported platform, skipping.",
                experiment_accession_code,
                survey_job=self.survey_job.id)
            return None, []

        return experiment, samples
