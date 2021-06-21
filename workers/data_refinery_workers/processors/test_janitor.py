import os
from unittest.mock import MagicMock, patch

from django.test import TestCase, tag

from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    Experiment,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    ProcessorJob,
    ProcessorJobDatasetAssociation,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleComputedFileAssociation,
    SampleResultAssociation,
    SurveyJob,
)
from data_refinery_common.utils import get_env_variable
from data_refinery_workers.processors import janitor

LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
JOBS = 10


def prepare_job():

    # Create 10 job directories
    for i in range(0, JOBS):

        os.makedirs(LOCAL_ROOT_DIR + "/processor_job_" + str(i), exist_ok=True)

        # These live on prod volumes at locations such as:
        # /var/ebs/SRP057116/SRR1972985/SRR1972985.sra
        os.makedirs(LOCAL_ROOT_DIR + "/SRP" + str(i), exist_ok=True)
        os.makedirs(LOCAL_ROOT_DIR + "/SRP" + str(i) + "/SRR" + str(i), exist_ok=True)

        sample = Sample()
        sample.accession_code = "SRR" + str(i)
        sample.save()

        cr = ComputationalResult()
        cr.save()

        cf = ComputedFile()
        cf.result = cr
        cf.size_in_bytes = 666
        cf.save()

        scfa = SampleComputedFileAssociation()
        scfa.sample = sample
        scfa.computed_file = cf
        scfa.save()

    # Create a job out of the range with index in it to make sure we
    # don't delete index directories since that's where transcriptome
    # indices get downloaded to.
    os.makedirs(LOCAL_ROOT_DIR + "/processor_job_" + str(JOBS + 1) + "_index", exist_ok=True)

    os.makedirs(LOCAL_ROOT_DIR + "/SRP" + str(JOBS + 1) + "/SRR" + str(JOBS + 1), exist_ok=True)
    sample = Sample()
    sample.accession_code = "SRR" + str(JOBS + 1)
    sample.save()

    # Save two jobs so that we trigger two special circumstances, one
    # where the job is still running and the other where querying
    # nomad raises an exception.
    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.nomad_job_id = "running_job"
    pj.save()

    pj = ProcessorJob()
    pj.pipeline_applied = "SALMON"
    pj.nomad_job_id = "missing_job"
    pj.save()

    pj = ProcessorJob()
    pj.pipeline_applied = "JANITOR"
    pj.save()

    return pj


class JanitorTestCase(TestCase):
    @tag("janitor")
    @patch("data_refinery_workers.processors.janitor.Nomad")
    def test_janitor(self, mock_nomad):
        """ Main tester. """

        def mock_get_job(job_id: str):
            if job_id == "running_job":
                return {"Status": "running"}
            else:
                return {"Status": "dead"}

        def mock_init_nomad(host, port=0, timeout=0):
            ret_value = MagicMock()
            ret_value.job = MagicMock()
            ret_value.job.get_job = MagicMock()
            ret_value.job.get_job.side_effect = mock_get_job
            return ret_value

        mock_nomad.side_effect = mock_init_nomad
        job = prepare_job()
        final_context = janitor.run_janitor(job.pk)

        for i in range(0, JOBS):
            # The job with id 1 should appear running.
            if i == 1:
                self.assertTrue(os.path.exists(LOCAL_ROOT_DIR + "/processor_job_" + str(i)))
            else:
                self.assertFalse(os.path.exists(LOCAL_ROOT_DIR + "/processor_job_" + str(i)))

            self.assertFalse(os.path.exists(LOCAL_ROOT_DIR + "/SRP" + str(i) + "/SRR" + str(i)))

        self.assertTrue(os.path.exists(LOCAL_ROOT_DIR + "/processor_job_11_index"))
        self.assertTrue(
            os.path.exists(LOCAL_ROOT_DIR + "/SRP" + str(JOBS + 1) + "/SRR" + str(JOBS + 1))
        )

        # Deleted all the working directories except for the one that's still running.
        self.assertEqual(len(final_context["deleted_items"]), (JOBS * 2) - 1)
