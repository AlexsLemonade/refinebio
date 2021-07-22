import datetime
from unittest.mock import patch

from django.test import TestCase

import vcr

from data_refinery_common.models import (
    DownloaderJob,
    Experiment,
    Organism,
    Sample,
    SurveyJob,
    SurveyJobKeyValue,
)
from data_refinery_foreman.surveyor.sra import SraSurveyor
from data_refinery_foreman.surveyor.surveyor import run_job

EXPERIMENT_ACCESSION = "DRX001563"
RUN_ACCESSION = "DRR002116"
SAMPLE_ACCESSION = "DRS001521"
STUDY_ACCESSION = "DRP000595"
SUBMISSION_ACCESSION = "DRA000567"


class SraSurveyorTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        self.survey_job = survey_job

        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value="DRR002116"
        )
        key_value_pair.save()

        # Insert the organism into the database so the model doesn't call the
        # taxonomy API to populate it.
        organism = Organism(name="HOMO_SAPIENS", taxonomy_id=9606, is_scientific_name=True)
        organism.save()

        organism1 = Organism(name="GALLUS_GALLUS", taxonomy_id=9031, is_scientific_name=True)
        organism1.save()

        organism2 = Organism(name="DANIO_RERIO", taxonomy_id=7955, is_scientific_name=True)
        organism2.save()

    def tearDown(self):
        DownloaderJob.objects.all().delete()
        SurveyJobKeyValue.objects.all().delete()
        SurveyJob.objects.all().delete()

    def test_survey(self):
        """A Simple test of the SRA surveyor.
        """
        sra_surveyor = SraSurveyor(self.survey_job)
        sra_surveyor.discover_experiment_and_samples()

        samples = Sample.objects.all()

        # We are expecting this to discover 1 sample.
        self.assertEqual(samples.count(), 1)
        # Confirm the sample's protocol_info
        experiment = Experiment.objects.all().first()
        self.assertEqual(
            samples.first().protocol_info[0]["Description"], experiment.protocol_description
        )

    @vcr.use_cassette("/home/user/data_store/cassettes/surveyor.sra.srp_survey.yaml")
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_srp_survey(self, mock_send_job):
        """A slightly harder test of the SRA surveyor.
        """
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value="SRP068364"
        )
        key_value_pair.save()

        sra_surveyor = SraSurveyor(survey_job)
        experiment, samples = sra_surveyor.discover_experiment_and_samples()
        self.assertEqual(experiment.accession_code, "SRP068364")
        self.assertEqual(experiment.alternate_accession_code, "GSE76780")
        self.assertEqual(len(samples), 4)

        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value="SRP111553"
        )
        key_value_pair.save()

        sra_surveyor = SraSurveyor(survey_job)
        experiment, samples = sra_surveyor.discover_experiment_and_samples()

        self.assertEqual(experiment.accession_code, "SRP111553")
        self.assertEqual(experiment.alternate_accession_code, "GSE101204")
        self.assertEqual(len(samples), 16)  # 8 samples with 2 runs each

        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value="DRP003977"
        )
        key_value_pair.save()

        sra_surveyor = SraSurveyor(survey_job)
        experiment, samples = sra_surveyor.discover_experiment_and_samples()

        self.assertEqual(experiment.accession_code, "DRP003977")
        self.assertEqual(experiment.alternate_accession_code, None)
        self.assertEqual(len(samples), 9)

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.sra.survey_file_report.yaml", record_mode="all"
    )
    @patch("data_refinery_foreman.surveyor.external_source.send_job")
    def test_survey_file_report(self, mock_send_job):
        """A slightly harder test of the SRA surveyor.
        """
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value="SRP048683"
        )
        key_value_pair.save()

        sra_surveyor = SraSurveyor(survey_job)
        experiment, samples = sra_surveyor.discover_experiment_and_samples()

        self.assertEqual(experiment.accession_code, "SRP048683")
        self.assertEqual(len(samples), 12)

        # Just check one file for one sample's expected file size/md5
        for sample in samples:
            if sample.accession_code == "SRR1603661":
                for original_file in sample.original_files.all():
                    if original_file.source_filename == "SRR1603661_1.fastq.gz":
                        self.assertEqual(
                            original_file.expected_md5, "502a9a482bfa5aa75865ccc0105ad13c"
                        )
                        self.assertEqual(original_file.expected_size_in_bytes, 6751980628)

    @vcr.use_cassette("/home/user/data_store/cassettes/surveyor.sra.survey_nonexistant.yaml")
    def test_nonexistant_srp_survey(self):
        """Try surveying an accession that does not exist
        """
        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value="ERP006216"
        )
        key_value_pair.save()

        run_job(survey_job)

        survey_job.refresh_from_db()
        self.assertFalse(survey_job.success)
        self.assertEqual(survey_job.failure_reason, "No experiment found.")

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.sra.arrayexpress_alternate_accession.yaml"
    )
    def test_arrayexpress_alternate_accession(self):
        """ Make sure that ENA experiments correctly detect their ArrayExpress alternate accession
        """

        survey_job = SurveyJob(source_type="SRA")
        survey_job.save()
        key_value_pair = SurveyJobKeyValue(
            survey_job=survey_job, key="experiment_accession_code", value="ERP108370"
        )
        key_value_pair.save()

        sra_surveyor = SraSurveyor(survey_job)
        experiment, _ = sra_surveyor.discover_experiment_and_samples()

        self.assertEqual(experiment.accession_code, "ERP108370")
        self.assertEqual(experiment.alternate_accession_code, "E-MTAB-6681")

    @vcr.use_cassette(
        "/home/user/data_store/cassettes/surveyor.sra.metadata_is_gathered_correctly.yaml"
    )
    def test_metadata_is_gathered_correctly(self):

        metadata = SraSurveyor.gather_all_metadata("DRR002116")

        self.assertEqual(metadata["broker_name"], "DDBJ")
        self.assertEqual(metadata["center_name"], "RIKEN_CDB")
        self.assertEqual(metadata["ena-base-count"], "158881910957")
        self.assertEqual(metadata["ena-spot-count"], "1371813555")
        self.assertEqual(metadata["experiment_accession"], "DRX001563")
        self.assertEqual(
            metadata["experiment_design_description"],
            ("Experiment for mRNAseq of chicken at stage " "HH16 (biological replicate 1)"),
        )
        self.assertEqual(
            metadata["experiment_title"],
            ("Illumina HiSeq 2000 sequencing; " "Exp_Gg_HH16_1_embryo_mRNAseq"),
        )
        self.assertEqual(
            metadata["lab_name"],
            (
                "Group for Morphological Evolution, Center for Developmental "
                "Biology, Kobe Institute, RIKEN"
            ),
        )
        self.assertEqual(metadata["library_layout"], "SINGLE")
        self.assertEqual(metadata["library_name"], "Gg_HH16_1_embryo_mRNAseq")
        self.assertEqual(metadata["library_selection"], "RANDOM")
        self.assertEqual(metadata["library_source"], "TRANSCRIPTOMIC")
        self.assertEqual(metadata["library_strategy"], "RNA-Seq")
        self.assertEqual(metadata["organism_id"], "9031")
        self.assertEqual(metadata["organism_name"], "GALLUS GALLUS")
        self.assertEqual(metadata["platform_instrument_model"], "Illumina HiSeq 2000")
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
        self.assertEqual(
            metadata["sample_sample_comment"],
            ("mRNAseq of chicken at stage HH16 (biological " "replicate 1)"),
        )
        self.assertEqual(metadata["sample_sample_name"], "DRS001521")
        self.assertEqual(metadata["sample_title"], "Gg_HH16_1_embryo_mRNAseq")
        self.assertEqual(metadata["spot_length"], "100")
        self.assertEqual(metadata["study_ena_first_public"], "2013-07-19")
        self.assertEqual(metadata["study_ena_last_update"], "2015-06-22")
        self.assertEqual(metadata["study_accession"], "DRP000595")
        self.assertEqual(metadata["submission_accession"], "DRA000567")
        self.assertEqual(
            metadata["submission_comment"],
            (
                "Time course gene expression profiles of turtle "
                "(Pelodiscus sinensis) and chicken (Gallus gallus) "
                "embryos were examined. Whole transcriptome of turtle "
                "was also determined by uding stranded sequencing "
                "methods."
            ),
        )
        self.assertEqual(metadata["submission_title"], "Submitted by RIKEN_CDB on 19-JUL-2013")

        ncbi_url = SraSurveyor._build_ncbi_file_url(metadata["run_accession"])
        self.assertTrue(
            ncbi_url
            in [
                "anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/DRR/DRR002/DRR002116/DRR002116.sra",
                "anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/DRR/DRR002/DRR002116/DRR002116.sra",
                "dbtest@sra-download.ncbi.nlm.nih.gov:data/sracloud/traces/dra0/DRR/000002/DRR002116",
            ]
        )

    def test_sra_metadata_is_harmonized(self):
        metadata = SraSurveyor.gather_all_metadata("SRR3098582")
        sample = Sample()
        SraSurveyor._apply_harmonized_metadata_to_sample(sample, metadata)
        self.assertEqual(sample.treatment, "biliatresone")
        self.assertEqual(sample.subject, "liver")
        self.assertEqual(sample.specimen_part, "liver")

        experiment = Experiment()
        SraSurveyor._apply_metadata_to_experiment(experiment, metadata)
        self.assertEqual(
            experiment.title,
            "Transcriptional profiling through RNA-seq of zebrafish larval"
            " liver after exposure to biliatresone, a biliary toxin.",
        )
        self.assertEqual(experiment.source_first_published, datetime.date(2017, 9, 25))
        self.assertEqual(experiment.source_last_modified, datetime.date(2017, 9, 25))
