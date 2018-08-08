"""
This command will remove experiments and all data that is uniquely associated with them.
(i.e. it's only associated with these experiment and not others.)
This can be used so that a surveyor job that went wrong can be rerun.
"""

from django.core.management.base import BaseCommand
from data_refinery_common.models import *
from data_refinery_common.logging import get_and_configure_logger

logger = get_and_configure_logger(__name__)

def purge_experiment(accession: str) -> None:
    experiment = Experiment.objects.filter(accession_code=accession)[0]
    ExperimentAnnotation.objects.filter(experiment=experiment).delete()
    experiment_sample_assocs = ExperimentSampleAssociation.objects.filter(experiment=experiment)
    samples = Sample.objects.filter(id__in=experiment_sample_assocs.values('sample_id'))

    uniquely_assoced_sample_ids = []
    deletable_experiment_sample_assocs = []
    for sample in samples:
        assoc_query = ExperimentSampleAssociation.objects.filter(sample=sample)
        if assoc_query.count() == 1:
            uniquely_assoced_sample_ids.append(sample.id)
            SampleAnnotation.objects.filter(sample=sample).delete()
            deletable_experiment_sample_assocs.append(assoc_query[0].id)

    comp_result_sample_assocs = SampleResultAssociation.objects.filter(sample_id__in=uniquely_assoced_sample_ids)
    comp_results = ComputationalResult.objects.filter(id__in=comp_result_sample_assocs.values('computational_result_id'))

    uniquely_assoced_comp_result_ids = []
    deletable_sample_result_assocs = []
    for comp_result in comp_results:
        extra_assocs = SampleResultAssociation.objects.filter(
            computational_result=comp_result,
            sample_id__not_in=uniquely_assoced_sample_ids)
        if extra_assocs.count == 0:
            uniquely_assoced_comp_result_ids.append(comp_result.id)

            unique_associations = SampleResultAssociation.objects.filter(result=comp_result)
            for assoc in unique_associations:
                deletable_sample_result_assocs.append(assoc)

            ComputationalResultAnnotation.objects.filter(result=comp_result).delete()
            computed_files = ComputedFile.objects.filter(result=comp_result)
            for computed_file in computed_files:
                computed_file.delete_local_file()
                computed_file.delete()

    og_file_sample_assocs = OriginalFileSampleAssociation.objects.filter(sample_id__in=uniquely_assoced_sample_ids)
    original_files = OriginalFile.objects.filter(id__in=og_file_sample_assocs.values('original_file_id'))

    uniquely_assoced_og_file_ids = []
    deletable_og_file_assocs = []
    for original_file in original_files:
        extra_assocs = OriginalFileSampleAssociation.objects.filter(
            original_file=original_file,
            sample_id__not_in=uniquely_assoced_sample_ids)
        if extra_assocs.count == 0:
            uniquely_assoced_og_file_ids.append(original_file.id)

            unique_associations = OriginalFilesampleAssociation.objects.filter(original_file=original_file)
            for assoc in unique_associations:
                deletable_og_file_assocs.append(assoc)

    og_file_dj_assocs = DownloaderJobOriginalFileAssociation.objects.filter(original_file_id__in=uniquely_assoced_og_file_ids)
    downloader_jobs = DownloaderJob.objects.filter(id__in=og_file_dj_assocs.values('downloader_job_id'))

    for downloader_job in downloader_jobs:
        extra_assocs = DownloaderJobOriginalFileAssociation.objects.filter(
            downloader_job=downloader_job,
            original_file_id__not_in=uniquely_assoced_og_file_ids)
        if extra_assocs.count == 0:
            DownloaderJobOriginalFileAssociation.objects.filter(downloader_job=downloader_job).delete()
            downloader_job.delete()

    for processor_job in processor_jobs:
        extra_assocs = ProcessorJobOriginalFileAssociation.objects.filter(
            processor_job=processor_job,
            original_file_id__not_in=uniquely_assoced_og_file_ids)
        if extra_assocs.count == 0:
            ProcessorJobOriginalFileAssociation.objects.filter(processor_job=processor_job).delete()
            processor_job.delete()

    for og_file_assoc in deletable_og_file_assocs:
        og_file_assoc.delete()

    for comp_result_assoc in deletable_sample_result_assocs:
        comp_result_assoc.delete()

    for sample_assoc in deletable_experiment_sample_assocs:
        sample_assoc.delete()

    OriginalFiles.objects.filter(id__in=uniquely_assoced_og_file_ids).delete()
    ComputationalResult.objects.filter(id__in=uniquely_assoced_comp_result_ids).delete()
    Samples.objects.filter(id__in=uniquely_assoced_sample_ids).delete()
    experiment.delete()




class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument(
            "--accession",
            help=("An experiment accession code to survey, download, and process."))
        parser.add_argument(
            "--file",
            type=str,
            help=("An optional file listing accession codes. s3:// URLs are also accepted.")
        )

    def handle(self, *args, **options):
        if options["accession"] is None and options['file'] is None:
            logger.error("You must specify an experiment accession or file.")
            sys.exit(1)

        accessions = []
        if options["file"]:
            if 's3://' in options["file"]:
                bucket, key = parse_s3_url(options["file"])
                s3 = boto3.resource('s3')
                try:
                    filepath = "/tmp/input_" + str(uuid.uuid4()) + ".txt"
                    s3.Bucket(bucket).download_file(key, filepath)
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        logger.error("The remote file does not exist.")
                        raise
                    else:
                        raise
            else:
                filepath = options["file"]

            with open(filepath) as file:
                for accession in file:
                    accessions.append(accession)
        else:
            accessions.append(options['accession'])

        for accession in accessions:
            logger.info("Purging Experiment with accession: %s", accession)
            try:
                purge_experiment(accession)
            except Exception as e:
                logger.exception("Exception caught while purging experiment with accession: %s", accession)

        logger.info("Purged %s experiments.", str(len(accessions)))
