import GEOparse
import dateutil.parser
import logging
import requests
import shutil

from re import sub, split, match
from typing import List, Dict

from data_refinery_common.job_lookup import ProcessorPipeline, Downloaders
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
    get_supported_rnaseq_platforms,
)
from data_refinery_foreman.surveyor import utils, harmony
from data_refinery_foreman.surveyor.external_source import ExternalSourceSurveyor


logger = get_and_configure_logger(__name__)
GEOparse.logger.set_verbosity("WARN")


UNKNOWN = "UNKNOWN"


class GeoSurveyor(ExternalSourceSurveyor):

    """Surveys NCBI GEO for data.

    Implements the ExternalSourceSurveyor interface.
    """

    def source_type(self):
        return Downloaders.GEO.value

    def get_temp_path(self):
        return '/tmp/' + str(self.survey_job.id) + '/'

    def set_platform_properties(self,
                                sample_object: Sample,
                                sample_metadata: Dict,
                                gse: GEOparse.GSM) -> Sample:
        """Sets platform-related properties on `sample_object`.

        Uses metadata from `gse` to populate platform_name,
        platform_accession_code, and technology on `sample_object`.
        """

        # Determine platform information
        external_accession = get_normalized_platform(gse.metadata.get('platform_id', [UNKNOWN])[0])

        if external_accession == UNKNOWN:
            sample_object.platform_accession_code = UNKNOWN
            sample_object.platform_name = UNKNOWN
            sample_object.manufacturer = UNKNOWN
            # If this sample is Affy, we potentially can extract the
            # platform information from the .CEL file. If it's not we
            # can't do anything. Therefore assume the technology is
            # microarray when we have no platform information.
            sample_object.technology = "MICROARRAY"
            return sample_object

        platform_accession_code = UNKNOWN

        gpl = GEOparse.get_GEO(external_accession, destdir=self.get_temp_path(), how="brief", silent=True)
        platform_title = gpl.metadata.get("title", [UNKNOWN])[0]

        # Check if this is a supported microarray platform.
        for platform in get_supported_microarray_platforms():
            if platform["external_accession"] == external_accession:
                platform_accession_code = platform["platform_accession"]

        if platform_accession_code != UNKNOWN:
            # It's a supported microarray platform.

            # We are using the brain array package as the platform accession code,
            # so, for instance, GPL3213 becomes 'chicken'.
            sample_object.platform_accession_code = platform_accession_code
            sample_object.technology = "MICROARRAY"
            try:

                # Related: https://github.com/AlexsLemonade/refinebio/issues/354
                # If it's Affy we can get a readable name:
                sample_object.platform_name = get_readable_affymetrix_names()[
                    platform_accession_code]
                sample_object.manufacturer = "AFFYMETRIX"

                # Sometimes Affymetrix samples have weird channel
                # protocol metadata, so if we find that it's
                # Affymetrix return it now. Example: GSE113945
                return sample_object
            except KeyError:
                # Otherwise we'll use what we've got.
                sample_object.platform_name = platform_title

            # Determine manufacturer

            # Sometimes this field is a list, other times it's not.
            # Example of it being a list: GSE113945
            channel1_temp = sample_metadata.get('label_protocol_ch1', "")
            if type(channel1_temp) != list:
                channel1_protocol = channel1_temp.upper()
            elif len(channel1_temp) > 0:
                channel1_protocol = channel1_temp[0].upper()
            else:
                channel1_protocol = ""

            if ('AGILENT' in channel1_protocol):
                sample_object.manufacturer = "AGILENT"
            elif ('ILLUMINA' in channel1_protocol):
                sample_object.manufacturer = "ILLUMINA"
            elif ('AFFYMETRIX' in channel1_protocol):
                sample_object.manufacturer = "AFFYMETRIX"
            else:
                sample_object.manufacturer = UNKNOWN

            return sample_object

        # Check to see if this is a supported RNASeq technology:

        # GEO RNASeq platform titles often have organisms appended to
        # an otherwise recognizable platform. The list of supported
        # RNASeq platforms isn't long, so see if any of them are
        # contained within what GEO gave us.
        # Example: GSE69572 has a platform title of:
        # 'Illumina Genome Analyzer IIx (Glycine max)'
        # Which should really just be 'Illumina Genome Analyzer IIx'
        # because RNASeq platforms are organism agnostic.  However,
        # the platforms 'Illumina Genome Analyzer' and 'Illumina
        # Genome Analyzer II' would also be matched, so make sure that
        # the longest platform names are tested first:
        sorted_platform_list = get_supported_rnaseq_platforms().copy()
        sorted_platform_list.sort(key=len, reverse=True)

        for platform in sorted_platform_list:
            if platform.upper() in platform_title.upper():
                sample_object.technology = "RNA-SEQ"
                sample_object.platform_name = platform
                # We just use RNASeq platform titles as accessions
                sample_object.platform_accession_code = platform

                if "ILLUMINA" in sample_object.platform_name.upper():
                    sample_object.manufacturer = "ILLUMINA"
                elif "NEXTSEQ" in sample_object.platform_name.upper():
                    sample_object.manufacturer = "NEXTSEQ"
                elif "ION TORRENT" in sample_object.platform_name.upper():
                    sample_object.manufacturer = "ION_TORRENT"
                else:
                    sample_object.manufacturer = UNKNOWN

                return sample_object

        # If we've made it this far, we don't know what this platform
        # is, therefore we can't know what its technology is. What we
        # do know is what GEO said was it's platform's accession and
        # title are, and that it's unsupported.
        sample_object.platform_name = platform_title
        sample_object.platform_accession_code = external_accession
        sample_object.technology = UNKNOWN
        sample_object.manufacturer = UNKNOWN

        return sample_object

    def get_miniml_url(self, experiment_accession_code):
        """ Build the URL for the MINiML files for this accession code.
        ex:
        'GSE68061' -> 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68061/miniml/GSE68061_family.xml.tgz'

        """
        geo = experiment_accession_code.upper()
        geotype = geo[:3]
        range_subdir = sub(r"\d{1,3}$", "nnn", geo)

        min_url_template = ("ftp://ftp.ncbi.nlm.nih.gov/geo/"
                            "series/{range_subdir}/{record}/miniml/{record_file}")
        min_url = min_url_template.format(range_subdir=range_subdir,
                                          record=geo,
                                          record_file="%s_family.xml.tgz" % geo)

        return min_url

    @staticmethod
    def get_sample_protocol_info(sample_metadata, sample_accession_code):
        protocol_info = dict()
        if 'extract_protocol_ch1' in sample_metadata:
            protocol_info['Extraction protocol'] = sample_metadata['extract_protocol_ch1']
        if 'label_protocol_ch1' in sample_metadata:
            protocol_info['Label protocol'] = sample_metadata['label_protocol_ch1']
        if 'hyb_protocol' in sample_metadata:
            protocol_info['Hybridization protocol'] = sample_metadata['hyb_protocol']
        if 'scan_protocol' in sample_metadata:
            protocol_info['Scan protocol'] = sample_metadata['scan_protocol']
        if 'data_processing' in sample_metadata:
            protocol_info['Data processing'] = sample_metadata['data_processing']

        protocol_info['Reference'] = ('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
                                              + sample_accession_code)
        return protocol_info

    def create_experiment_and_samples_from_api(self, experiment_accession_code) -> (Experiment, List[Sample]):
        """ The main surveyor - find the Experiment and Samples from NCBI GEO.

        Uses the GEOParse library, for which docs can be found here: https://geoparse.readthedocs.io/en/latest/usage.html#working-with-geo-objects

        """
        # Cleaning up is tracked here: https://github.com/guma44/GEOparse/issues/41
        gse = GEOparse.get_GEO(experiment_accession_code, destdir=self.get_temp_path(), how="brief", silent=True)
        preprocessed_samples = harmony.preprocess_geo(gse.gsms.items())
        harmonized_samples = harmony.harmonize(preprocessed_samples)

        # Create the experiment object
        try:
            experiment_object = Experiment.objects.get(accession_code=experiment_accession_code)
            logger.debug("Experiment %s already exists, skipping object creation.",
                         experiment_accession_code,
                         survey_job=self.survey_job.id)
        except Experiment.DoesNotExist:
            experiment_object = Experiment()
            experiment_object.accession_code = experiment_accession_code
            experiment_object.source_url = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
                                            + experiment_accession_code)
            experiment_object.source_database = "GEO"
            experiment_object.title = gse.metadata.get('title', [''])[0]
            experiment_object.description = gse.metadata.get('summary', [''])[0]

            # Source doesn't provide time information, assume midnight.
            submission_date = gse.metadata["submission_date"][0] + " 00:00:00 UTC"
            experiment_object.source_first_published = dateutil.parser.parse(submission_date)
            last_updated_date = gse.metadata["last_update_date"][0] + " 00:00:00 UTC"
            experiment_object.source_last_updated = dateutil.parser.parse(last_updated_date)

            unique_institutions = list(set(gse.metadata["contact_institute"]))
            experiment_object.submitter_institution = ", ".join(unique_institutions)
            experiment_object.pubmed_id = gse.metadata.get("pubmed_id", [""])[0]

            # Scrape publication title and authorship from Pubmed
            if experiment_object.pubmed_id:
                pubmed_metadata = utils.get_title_and_authors_for_pubmed_id(experiment_object.pubmed_id)
                experiment_object.publication_title = pubmed_metadata[0]
                experiment_object.publication_authors = pubmed_metadata[1]

            experiment_object.save()

            experiment_annotation = ExperimentAnnotation()
            experiment_annotation.data = gse.metadata
            experiment_annotation.experiment = experiment_object
            experiment_annotation.is_ccdl = False
            experiment_annotation.save()

        # Okay, here's the situation!
        # Sometimes, samples have a direct single representation for themselves.
        # Othertimes, there is a single file with references to every sample in it.
        created_samples = []
        for sample_accession_code, sample in gse.gsms.items():

            try:
                sample_object = Sample.objects.get(accession_code=sample_accession_code)
                logger.debug(
                    "Sample %s from experiment %s already exists, skipping object creation.",
                         sample_accession_code,
                         experiment_object.accession_code,
                         survey_job=self.survey_job.id)

                # Associate it with the experiment, but since it
                # already exists it already has original files
                # associated with it and it's already been downloaded,
                # so don't add it to created_samples.
                ExperimentSampleAssociation.objects.get_or_create(
                    experiment=experiment_object, sample=sample_object)
            except Sample.DoesNotExist:
                organism = Organism.get_object_for_name(sample.metadata['organism_ch1'][0].upper())

                sample_object = Sample()
                sample_object.source_database = "GEO"
                sample_object.accession_code = sample_accession_code
                sample_object.organism = organism

                # If data processing step, it isn't raw.
                sample_object.has_raw = not sample.metadata.get('data_processing', None)

                ExperimentOrganismAssociation.objects.get_or_create(
                    experiment=experiment_object, organism=organism)
                sample_object.title = sample.metadata['title'][0]

                self.set_platform_properties(sample_object, sample.metadata, gse)

                # Directly assign the harmonized properties
                harmonized_sample = harmonized_samples[sample_object.title]
                for key, value in harmonized_sample.items():
                    setattr(sample_object, key, value)

                # Sample-level protocol_info
                sample_object.protocol_info = self.get_sample_protocol_info(sample.metadata,
                                                                            sample_accession_code)

                sample_object.save()
                created_samples.append(sample_object)
                logger.debug("Created Sample: " + str(sample_object))

                # Now that we've determined the technology at the
                # sample level, we can set it at the experiment level,
                # just gotta make sure to only do it once.
                # XXX: Do we want to remove this field? It's not going
                # to be accurate in cases where an experiment has both
                # technologies.
                if not experiment_object.technology:
                    experiment_object.technology = sample_object.technology
                    experiment_object.save()

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
                    if "idat.gz" in supplementary_file_url.lower():
                        continue
                    if "chp.gz" in supplementary_file_url.lower():
                        continue
                    if "ndf.gz" in supplementary_file_url.lower():
                        continue
                    if "pos.gz" in supplementary_file_url.lower():
                        continue
                    if "pair.gz" in supplementary_file_url.lower():
                        continue
                    if "gff.gz" in supplementary_file_url.lower():
                        continue

                    # Sometimes, we are lied to about the data processing step.
                    lower_file_url = supplementary_file_url.lower()
                    if '.cel' in lower_file_url \
                    or ('_non_normalized.txt' in lower_file_url) \
                    or ('_non-normalized.txt' in lower_file_url) \
                    or ('-non-normalized.txt' in lower_file_url) \
                    or ('-non_normalized.txt' in lower_file_url):
                        sample_object.has_raw = True
                        sample_object.save()

                    original_file = OriginalFile.objects.get_or_create(
                            source_url = supplementary_file_url,
                            source_filename = supplementary_file_url.split('/')[-1],
                            has_raw = sample_object.has_raw,
                            is_archive = True
                        )[0]

                    original_file_sample_association = OriginalFileSampleAssociation.objects.get_or_create(
                            original_file = original_file,
                            sample = sample_object
                        )

                ExperimentSampleAssociation.objects.get_or_create(
                    experiment=experiment_object, sample=sample_object)

        # These supplementary files _may-or-may-not_ contain the type of raw data we can process.
        for experiment_supplement_url in gse.metadata.get('supplementary_file', []):

            original_file = OriginalFile.objects.get_or_create(
                    source_url = experiment_supplement_url,
                    source_filename = experiment_supplement_url.split('/')[-1],
                    has_raw = sample_object.has_raw,
                    is_archive = True
                )[0]

            logger.debug("Created OriginalFile: " + str(original_file))

            lower_supplement_url = experiment_supplement_url.lower()
            if ('_non_normalized.txt' in lower_supplement_url) \
            or ('_non-normalized.txt' in lower_supplement_url) \
            or ('-non-normalized.txt' in lower_supplement_url) \
            or ('-non_normalized.txt' in lower_supplement_url):
                for sample_object in created_samples:
                    sample_object.has_raw = True
                    sample_object.save()

                    OriginalFileSampleAssociation.objects.get_or_create(
                        sample=sample_object, original_file=original_file)

            # Delete this Original file if it isn't being used.
            if OriginalFileSampleAssociation.objects.filter(original_file=original_file).count() == 0:
                original_file.delete()

        # These are the Miniml/Soft/Matrix URLs that are always(?) provided.
        # GEO describes different types of data formatting as "families"
        family_url = self.get_miniml_url(experiment_accession_code)
        miniml_original_file = OriginalFile.objects.get_or_create(
                source_url = family_url,
                source_filename = family_url.split('/')[-1],
                has_raw = sample_object.has_raw,
                is_archive = True
            )[0]
        for sample_object in created_samples:
            # We don't need a .txt if we have a .CEL
            if sample_object.has_raw:
                continue
            OriginalFileSampleAssociation.objects.get_or_create(
                sample=sample_object, original_file=miniml_original_file)

        # Delete this Original file if it isn't being used.
        if OriginalFileSampleAssociation.objects.filter(original_file=miniml_original_file).count() == 0:
            miniml_original_file.delete()

        # Trash the temp path
        try:
            shutil.rmtree(self.get_temp_path())
        except Exception:
            # There was a problem during surveying so this didn't get created.
            # It's not a big deal.
            pass

        return experiment_object, created_samples

    def discover_experiment_and_samples(self) -> (Experiment, List[Sample]):
        """ Dispatches the surveyor, returns the results """

        experiment_accession_code = (
            SurveyJobKeyValue
            .objects
            .get(survey_job_id=self.survey_job.id,
                 key__exact="experiment_accession_code")
            .value
        )

        logger.debug("Surveying experiment with accession code: %s.",
                    experiment_accession_code,
                    survey_job=self.survey_job.id)

        experiment, samples = self.create_experiment_and_samples_from_api(
            experiment_accession_code)

        return experiment, samples
