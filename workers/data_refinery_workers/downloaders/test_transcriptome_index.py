import copy
import shutil
from typing import List
from unittest.mock import call, patch

from django.test import TestCase, tag

from data_refinery_common.job_lookup import ProcessorPipeline
from data_refinery_common.models import (
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    OriginalFile,
    ProcessorJob,
    SurveyJob,
)
from data_refinery_workers.downloaders import transcriptome_index


class DownloadTranscriptomeIndexTestCase(TestCase):
    def setUp(self):
        survey_job = SurveyJob(source_type="TRANSCRIPTOME_INDEX")
        survey_job.save()
        self.survey_job = survey_job
        self.gtf_download_url = (
            "ftp://ftp.ensemblgenomes.org/pub/release-37/plants/gtf/"
            "aegilops_tauschii/Aegilops_tauschii.ASM34733v1.37.gtf.gz"
        )
        self.fasta_download_url = (
            "ftp://ftp.ensemblgenomes.org/pub/release-37/plants/fasta/"
            "aegilops_tauschii/dna/Aegilops_tauschii.ASM34733v1.dna.toplevel.fa.gz"
        )

    @tag("downloaders")
    @patch("data_refinery_workers.downloaders.transcriptome_index.send_job")
    def test_download_file(self, mock_send_job):
        mock_send_job.return_value = None
        dlj = DownloaderJob()
        dlj.save()
        og = OriginalFile()
        og.source_filename = "Aegilops_tauschii.ASM34733v1.37.gtf.gz"
        og.source_url = self.gtf_download_url
        og.is_archive = True
        og.save()

        assoc = DownloaderJobOriginalFileAssociation()
        assoc.downloader_job = dlj
        assoc.original_file = og
        assoc.save()

        transcriptome_index.download_transcriptome(dlj.pk)
