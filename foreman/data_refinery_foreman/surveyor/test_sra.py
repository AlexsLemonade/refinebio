from unittest.mock import Mock, patch, call
from django.test import TestCase
from data_refinery_foreman.surveyor.sra import (
    SraSurveyor,
    ENA_METADATA_URL_TEMPLATE,
)
from data_refinery_foreman.surveyor.test_sra_xml import (
    EXPERIMENT_XML,
    RUN_XML,
    SAMPLE_XML,
    STUDY_XML,
    SUBMISSION_XML
)
from data_refinery_common.models import (
    SurveyJob,
    SurveyJobKeyValue
)

EXPERIMENT_ACCESSION = "DRX001563"
RUN_ACCESSION = "DRR002116"
SAMPLE_ACCESSION = "DRS001521"
STUDY_ACCESSION = "DRP000595"
SUBMISSION_ACCESSION = "DRA000567"


def mocked_requests_get(url):
    mock = Mock(ok=True)
    if url == ENA_METADATA_URL_TEMPLATE.format(EXPERIMENT_ACCESSION):
        mock.text = EXPERIMENT_XML
    elif url == ENA_METADATA_URL_TEMPLATE.format(RUN_ACCESSION):
        mock.text = RUN_XML
    elif url == ENA_METADATA_URL_TEMPLATE.format(SAMPLE_ACCESSION):
        mock.text = SAMPLE_XML
    elif url == ENA_METADATA_URL_TEMPLATE.format(STUDY_ACCESSION):
        mock.text = STUDY_XML
    elif url == ENA_METADATA_URL_TEMPLATE.format(SUBMISSION_ACCESSION):
        mock.text = SUBMISSION_XML
    else:
        raise Exception("Was not expecting the url: " + url)

    return mock


class SraSurveyorTestCase(TestCase):
    def test_get_next_accession(self):
        self.assertEqual(SraSurveyor.get_next_accession("DRR123456"), "DRR123457")
        self.assertEqual(SraSurveyor.get_next_accession("DRR1234567"), "DRR1234568")
        self.assertEqual(SraSurveyor.get_next_accession("DRR12345678"), "DRR12345679")
        self.assertEqual(SraSurveyor.get_next_accession("DRR123456789"), "DRR123456790")

    @patch.object(SraSurveyor, "_generate_batch")
    def test_discover_batches(self, mock_generate_batch):
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()

        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="start_accession",
                                           value="DRR012345")
        key_value_pair.save()
        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="end_accession",
                                           value="DRR012348")
        key_value_pair.save()

        sra_surveyor = SraSurveyor(survey_job)
        sra_surveyor.discover_batches()

        mock_generate_batch.assert_has_calls([
            call("DRR012345"),
            call("DRR012346"),
            call("DRR012347"),
            call("DRR012348")
        ])

    @patch('data_refinery_foreman.surveyor.sra.requests.get')
    def test_metadata_is_gathered_correctly(self, mock_get):
        mock_get.side_effect = mocked_requests_get

        metadata = SraSurveyor.gather_all_metadata("DRR002116")

        self.assertEqual(metadata["broker_name"], "DDBJ")
        self.assertEqual(metadata["center_name"], "RIKEN_CDB")
        self.assertEqual(metadata["ena-base-count"], "158881910957")
        self.assertEqual(metadata["ena-spot-count"], "1371813555")
        self.assertEqual(metadata["experiment_accession"], "DRX001563")
        self.assertEqual(metadata["experiment_design_description"],
                         ("Experiment for mRNAseq of chicken at stage "
                          "HH16 (biological replicate 1)"))
        self.assertEqual(metadata["experiment_title"],
                         ("Illumina HiSeq 2000 sequencing; "
                          "Exp_Gg_HH16_1_embryo_mRNAseq"))
        self.assertEqual(metadata["lab_name"],
                         ("Group for Morphological Evolution, Center for Developmental "
                          "Biology, Kobe Institute, RIKEN"))
        self.assertEqual(metadata["library_layout"], "SINGLE")
        self.assertEqual(metadata["library_name"], "Gg_HH16_1_embryo_mRNAseq")
        self.assertEqual(metadata["library_selection"], "RANDOM")
        self.assertEqual(metadata["library_source"], "TRANSCRIPTOMIC")
        self.assertEqual(metadata["library_strategy"], "RNA-Seq")
        self.assertEqual(metadata["organism_id"], "9031")
        self.assertEqual(metadata["organism_name"], "GALLUS GALLUS")
        self.assertEqual(metadata["platform_instrument_model"], "IlluminaHiSeq2000")
        self.assertEqual(metadata["read_spec_0_base_coord"], "1")
        self.assertEqual(metadata["read_spec_0_class"], "Application Read")
        self.assertEqual(metadata["read_spec_0_index"], "0")
        self.assertEqual(metadata["read_spec_0_type"], "Forward")
        self.assertEqual(metadata["run_accession"], "DRR002116")
        self.assertEqual(metadata["run_center"], "RIKEN_CDB")
        self.assertEqual(metadata["run_date"], "2011-09-01T00:00:00+09:00")
        self.assertEqual(metadata["run_ena_base_count"], "3256836000")
        self.assertEqual(metadata["run_ena_first_public"], "2013-07-19")
        self.assertEqual(metadata["run_ena_last_update"], "2017-08-11")
        self.assertEqual(metadata["run_ena_spot_count"], "32568360")
        self.assertEqual(metadata["sample_accession"], "DRS001521")
        self.assertEqual(metadata["sample_center_name"], "BioSample")
        self.assertEqual(metadata["sample_ena_base_count"], "3256836000")
        self.assertEqual(metadata["sample_ena_first_public"], "2013-07-20")
        self.assertEqual(metadata["sample_ena_last_update"], "2015-08-24")
        self.assertEqual(metadata["sample_ena_spot_count"], "32568360")
        self.assertEqual(metadata["sample_sample_comment"],
                         ("mRNAseq of chicken at stage HH16 (biological "
                          "replicate 1)"))
        self.assertEqual(metadata["sample_sample_name"], "DRS001521")
        self.assertEqual(metadata["sample_title"], "Gg_HH16_1_embryo_mRNAseq")
        self.assertEqual(metadata["spot_length"], "100")
        self.assertEqual(metadata["study_accession"], "DRP000595")
        self.assertEqual(metadata["submission_accession"], "DRA000567")
        self.assertEqual(metadata["submission_comment"],
                         ("Time course gene expression profiles of turtle "
                          "(Pelodiscus sinensis) and chicken (Gallus gallus) "
                          "embryos were examined. Whole transcriptome of turtle "
                          "was also determined by uding stranded sequencing "
                          "methods."))
        self.assertEqual(metadata["submission_title"], "Submitted by RIKEN_CDB on 19-JUL-2013")

    @patch('data_refinery_foreman.surveyor.sra.requests.get')
    def test_batch_created(self, mock_get):
        mock_get.side_effect = mocked_requests_get

        # Use same run accession for the start and end of the range to
        # achieve a length of 1
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="start_accession",
                                           value=RUN_ACCESSION)
        key_value_pair.save()
        key_value_pair = SurveyJobKeyValue(survey_job=survey_job,
                                           key="end_accession",
                                           value=RUN_ACCESSION)
        key_value_pair.save()

        surveyor = SraSurveyor(survey_job)

        self.assertTrue(surveyor.discover_batches())
        # With only a single run accession there should only be a
        # single batch.
        self.assertEqual(len(surveyor.batches), 1)

        batch = surveyor.batches[0]
        self.assertEqual(batch.survey_job.id, survey_job.id)
        self.assertEqual(batch.source_type, "SRA")
        self.assertEqual(batch.pipeline_required, "SALMON")
        self.assertEqual(batch.platform_accession_code, "IlluminaHiSeq2000")
        self.assertEqual(batch.experiment_accession_code, "DRX001563")
        self.assertEqual(batch.experiment_title, ("Illumina HiSeq 2000 sequencing; "
                                                  "Exp_Gg_HH16_1_embryo_mRNAseq"))
        self.assertEqual(batch.status, "NEW")
        self.assertEqual(batch.release_date, "2013-07-19")
        self.assertEqual(batch.last_uploaded_date, "2017-08-11")
        self.assertEqual(batch.organism_id, 9031)
        self.assertEqual(batch.organism_name, "GALLUS GALLUS")

        file = batch.files[0]
        self.assertEqual(file.size_in_bytes, -1)
        self.assertEqual(file.download_url, "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR002/DRR002116/DRR002116.fastq.gz")  # noqa
        self.assertEqual(file.raw_format, "fastq.gz")
        self.assertEqual(file.processed_format, "tar.gz")
        self.assertEqual(file.name, "DRR002116.fastq.gz")
        self.assertEqual(file.internal_location, "IlluminaHiSeq2000/SALMON")
