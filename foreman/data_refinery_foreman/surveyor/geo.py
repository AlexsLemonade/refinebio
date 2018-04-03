import GEOparse
import requests

from re import sub, split, match

from data_refinery_foreman.surveyor import utils
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

class GeoUnsupportedPlatformException(BaseException):
    pass

class GeoSurveyor(ExternalSourceSurveyor):
    """Surveys NCBI GEO for data.

    Implements the ExternalSourceSurveyor interface.
    """

    def get_raw_url(self, experiment_accession_code):
        """ """
        geo = experiment_accession_code.upper()
        geotype = geo[:3]
        range_subdir = sub(r"\d{1,3}$", "nnn", geo)

        # miniml_url_template = ("ftp://ftp.ncbi.nlm.nih.gov/geo/"
        #           "{root}/{range_subdir}/{record}/miniml/{record_file}")
        # miniml_url = miniml_url_template.format(root="series",
        #                     range_subdir=range_subdir,
        #                     record=geo,
        #                     record_file="%s_family.xml.tgz" % geo)

        raw_url_template = ("ftp://ftp.ncbi.nlm.nih.gov/geo/"
                  "{root}/{range_subdir}/{record}/suppl/{record_file}")
        raw_url = raw_url_template.format(root="series",
                            range_subdir=range_subdir,
                            record=geo,
                            record_file="%s_RAW.tar" % geo)

        return raw_url

    def create_experiment_and_samples_from_api(self, experiment_accession_code) -> (Experiment, List[Sample]):
        """ """

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
            experiment_object.name = gse.metadata['title'][0]
            experiment_object.description = gse.metadata['summary'][0]
            experiment_object.platform_name = gse.metadata["platform_id"][0] # TODO: Lookup GEO-GPL
            experiment_object.platform_accession_code = gse.metadata["platform_id"][0]
            experiment_object.source_first_published = gse.metadata["submission_date"][0]
            experiment_object.source_last_updated = gse.metadata["last_update_date"][0]
            experiment_object.submitter_institution = ", ".join(list(set(gse.metadata["contact_institute"])))
            experiment_object.pubmed_id = gse.metadata.get("pubmed_id", [None])[0]
            experiment_object.save()

            experiment_annotation = ExperimentAnnotation()
            experiment_annotation.data = gse.metadata
            experiment_annotation.experiment = experiment_object
            experiment_annotation.is_ccdl = False
            experiment_annotation.save()

            has_raw = True

        created_samples = []
        for sample_accession_code, sample in gse.gsms.items():

            try:
                sample_object = Sample.objects.get(accession_code=sample_accession_code)
                logger.error("Sample %s from experiment %s already exists, skipping object creation.",
                         sample_accession_code,
                         experiment_object.accession_code,
                         survey_job=self.survey_job.id)
                continue
            except Sample.DoesNotExist:
                organism = Organism.get_object_for_name(sample.metadata['organism_ch1'][0].upper())

                sample_object = Sample()
                sample_object.accession_code = sample_accession_code
                sample_object.organism = organism
                sample_object.save()

                sample_annotation = SampleAnnotation()
                sample_annotation.sample = sample
                sample_annotation.data = sample.metadata
                sample_annotation.is_ccdl = False
                sample_annotation.save()

                original_file = OriginalFile()
                original_file.sample = sample_object
                original_file.source_filename = sample_accession_code + '.txt.gz'
                original_file.source_url = get_raw_url(experiment_accession_code)
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

        return experiment, created_samples

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