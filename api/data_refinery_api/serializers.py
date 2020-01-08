from collections import defaultdict

from django.db.models import Count, Q
from rest_framework import serializers

import boto3
from django_elasticsearch_dsl_drf.serializers import DocumentSerializer

from data_refinery_common.models import (
    APIToken,
    CompendiumResult,
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    Dataset,
    DownloaderJob,
    Experiment,
    ExperimentAnnotation,
    Organism,
    OrganismIndex,
    OriginalFile,
    Processor,
    ProcessorJob,
    Sample,
    SampleAnnotation,
    SurveyJob,
)
from data_refinery_common.models.documents import ExperimentDocument

s3 = boto3.client("s3")


##
# Organism
##


class OrganismSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = (
            "name",
            "taxonomy_id",
        )


##
# Processor
##


class ProcessorSerializer(serializers.ModelSerializer):
    class Meta:
        model = Processor
        fields = ("id", "name", "version", "docker_image", "environment")


##
# Transcriptome Index
##


class OrganismIndexSerializer(serializers.ModelSerializer):

    organism_name = serializers.StringRelatedField(source="organism", read_only=True)
    download_url = serializers.SerializerMethodField()

    class Meta:
        model = OrganismIndex

        fields = (
            "id",
            "assembly_name",
            "organism_name",
            "source_version",
            "index_type",
            "salmon_version",
            "download_url",
            "result_id",
            "last_modified",
        )
        read_only_fields = fields

    def get_download_url(self, obj):
        computed_file = obj.get_computed_file()
        if computed_file is not None:
            return computed_file.s3_url
        return None


##
# Results
##


class DetailedExperimentSampleSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sample
        fields = (
            "accession_code",
            "platform_name",
            "pretty_platform",
            "technology",
            "is_processed",
        )


class ComputationalResultAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = ComputationalResultAnnotation
        fields = ("id", "data", "is_ccdl", "created_at", "last_modified")


class ComputedFileSerializer(serializers.ModelSerializer):
    class Meta:
        model = ComputedFile
        fields = (
            "id",
            "filename",
            "size_in_bytes",
            "is_smashable",
            "is_qc",
            "sha1",
            "s3_bucket",
            "s3_key",
            "created_at",
            "last_modified",
        )


class ComputedFileWithUrlSerializer(serializers.ModelSerializer):
    class Meta:
        model = ComputedFile
        fields = (
            "id",
            "filename",
            "size_in_bytes",
            "is_smashable",
            "is_qc",
            "sha1",
            "s3_bucket",
            "s3_key",
            "download_url",
            "created_at",
            "last_modified",
        )


class ComputationalResultSerializer(serializers.ModelSerializer):
    annotations = ComputationalResultAnnotationSerializer(
        many=True, source="computationalresultannotation_set"
    )
    processor = ProcessorSerializer(many=False)
    organism_index = OrganismIndexSerializer(many=False)
    files = ComputedFileSerializer(many=True, source="computedfile_set")

    class Meta:
        model = ComputationalResult
        fields = (
            "id",
            "commands",
            "processor",
            "is_ccdl",
            "annotations",
            "files",
            "organism_index",
            "time_start",
            "time_end",
            "created_at",
            "last_modified",
        )


class ComputationalResultWithUrlSerializer(ComputationalResultSerializer):
    files = ComputedFileWithUrlSerializer(many=True, source="computedfile_set")


class ComputationalResultNoFilesSerializer(serializers.ModelSerializer):
    annotations = ComputationalResultAnnotationSerializer(
        many=True, source="computationalresultannotation_set"
    )
    processor = ProcessorSerializer(many=False)
    organism_index = OrganismIndexSerializer(many=False)

    class Meta:
        model = ComputationalResult
        fields = (
            "id",
            "commands",
            "processor",
            "is_ccdl",
            "annotations",
            "organism_index",
            "time_start",
            "time_end",
            "created_at",
            "last_modified",
        )


class QNTargetSerializer(serializers.ModelSerializer):
    result = ComputationalResultNoFilesSerializer(many=False)

    class Meta:
        model = ComputedFile
        fields = (
            "id",
            "filename",
            "size_in_bytes",
            "is_qn_target",
            "sha1",
            "s3_bucket",
            "s3_key",
            "s3_url",
            "created_at",
            "last_modified",
            "result",
        )


class ComputedFileListSerializer(serializers.ModelSerializer):
    result = ComputationalResultNoFilesSerializer(many=False)
    samples = DetailedExperimentSampleSerializer(many=True)
    compendia_organism_name = serializers.CharField(
        source="compendia_organism__name", read_only=True
    )

    def __init__(self, *args, **kwargs):
        super(ComputedFileListSerializer, self).__init__(*args, **kwargs)
        if "context" in kwargs:
            # only include the field `download_url` if a valid token is specified
            # the token lookup happens in the view.
            if "token" not in kwargs["context"]:
                self.fields.pop("download_url")

    class Meta:
        model = ComputedFile
        fields = (
            "id",
            "filename",
            "samples",
            "size_in_bytes",
            "is_qn_target",
            "is_smashable",
            "is_qc",
            "is_compendia",
            "quant_sf_only",
            "compendia_version",
            "compendia_organism_name",
            "sha1",
            "s3_bucket",
            "s3_key",
            "s3_url",
            "download_url",
            "created_at",
            "last_modified",
            "result",
        )
        extra_kwargs = {
            "download_url": {
                "help_text": "This will contain an url to download the file. You must send a valid [token](#tag/token) in order to receive this."
            }
        }


class OriginalFileListSerializer(serializers.ModelSerializer):
    class Meta:
        model = OriginalFile
        fields = (
            "id",
            "filename",
            "samples",
            "size_in_bytes",
            "sha1",
            "samples",
            "processor_jobs",
            "downloader_jobs",
            "source_url",
            "is_archive",
            "source_filename",
            "has_raw",
            "created_at",
            "last_modified",
        )


##
# Samples
##


class SampleSerializer(serializers.ModelSerializer):
    organism = OrganismSerializer(many=False)

    class Meta:
        model = Sample
        fields = (
            "id",
            "title",
            "accession_code",
            "source_database",
            "organism",
            "platform_accession_code",
            "platform_name",
            "pretty_platform",
            "technology",
            "manufacturer",
            "protocol_info",
            "is_processed",
            "created_at",
            "last_modified",
        )


class SampleAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = SampleAnnotation
        fields = (
            "data",
            "is_ccdl",
            "created_at",
            "last_modified",
        )


class DetailedSamplesComputationalResultSerializer(serializers.ModelSerializer):
    processor = ProcessorSerializer(many=False)
    organism_index = OrganismIndexSerializer(many=False)

    class Meta:
        model = ComputationalResult
        fields = (
            "processor",
            "organism_index",
        )


class DetailedSampleSerializer(serializers.ModelSerializer):
    annotations = SampleAnnotationSerializer(many=True, source="sampleannotation_set")
    organism = OrganismSerializer(many=False)
    results = DetailedSamplesComputationalResultSerializer(many=True)

    class Meta:
        model = Sample
        fields = (
            "id",
            "title",
            "accession_code",
            "source_database",
            "organism",
            "platform_accession_code",
            "platform_name",
            "pretty_platform",
            "technology",
            "manufacturer",
            "protocol_info",
            "annotations",
            "results",
            "source_archive_url",
            "has_raw",
            "sex",
            "age",
            "specimen_part",
            "genotype",
            "disease",
            "disease_stage",
            "cell_line",
            "treatment",
            "race",
            "subject",
            "compound",
            "time",
            "is_processed",
            "created_at",
            "last_modified",
            "original_files",
            "computed_files",
        )


##
# Experiments
##


class ExperimentSerializer(serializers.ModelSerializer):
    organisms = serializers.StringRelatedField(many=True, read_only=True)
    platforms = serializers.ReadOnlyField()
    processed_samples = serializers.StringRelatedField(many=True)
    total_samples_count = serializers.IntegerField(read_only=True)
    sample_metadata = serializers.ReadOnlyField(source="get_sample_metadata_fields")
    technologies = serializers.ReadOnlyField(source="get_sample_technologies")
    pretty_platforms = serializers.ReadOnlyField()

    class Meta:
        model = Experiment
        fields = (
            "id",
            "title",
            "description",
            "accession_code",
            "alternate_accession_code",
            "source_database",
            "source_url",
            "platforms",
            "pretty_platforms",
            "processed_samples",
            "has_publication",
            "publication_title",
            "publication_doi",
            "publication_authors",
            "pubmed_id",
            "total_samples_count",
            "organisms",
            "submitter_institution",
            "created_at",
            "last_modified",
            "source_first_published",
            "source_last_modified",
            "sample_metadata",
            "technologies",
        )

    @staticmethod
    def setup_eager_loading(queryset):
        """ Perform necessary eager loading of data. """
        queryset = queryset.prefetch_related("samples").prefetch_related("organisms")

        # Multiple count annotations
        queryset = queryset.annotate(
            total_samples_count=Count("samples", unique=True),
            processed_samples_count=Count("samples", filter=Q(samples__is_processed=True)),
        )

        return queryset


class ExperimentAnnotationSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExperimentAnnotation
        fields = (
            "data",
            "is_ccdl",
            "created_at",
            "last_modified",
        )


class DetailedExperimentSerializer(serializers.ModelSerializer):
    annotations = ExperimentAnnotationSerializer(many=True, source="experimentannotation_set")
    samples = DetailedExperimentSampleSerializer(many=True)
    sample_metadata = serializers.ReadOnlyField(source="sample_metadata_fields")
    organisms = serializers.StringRelatedField(many=True, read_only=True)

    class Meta:
        model = Experiment
        fields = (
            "id",
            "title",
            "description",
            "annotations",
            "samples",
            "protocol_description",
            "accession_code",
            "source_database",
            "source_url",
            "has_publication",
            "publication_title",
            "publication_doi",
            "publication_authors",
            "pubmed_id",
            "source_first_published",
            "source_last_modified",
            "submitter_institution",
            "last_modified",
            "created_at",
            "organisms",
            "sample_metadata",
            "num_total_samples",
            "num_processed_samples",
            "num_downloadable_samples",
        )


class PlatformSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sample
        fields = (
            "platform_accession_code",
            "platform_name",
        )


class InstitutionSerializer(serializers.ModelSerializer):
    class Meta:
        model = Experiment
        fields = ("submitter_institution",)


##
# Files
##


class OriginalFileSerializer(serializers.ModelSerializer):
    class Meta:
        model = OriginalFile
        fields = (
            "id",
            "filename",
            "size_in_bytes",
            "sha1",
            "source_url",
            "source_filename",
            "is_downloaded",
            "is_archive",
            "has_raw",
            "created_at",
            "last_modified",
        )


##
# Jobs
##


class SurveyJobSerializer(serializers.ModelSerializer):
    class Meta:
        model = SurveyJob
        fields = (
            "id",
            "source_type",
            "success",
            "start_time",
            "end_time",
            "created_at",
            "last_modified",
        )


class DownloaderJobSerializer(serializers.ModelSerializer):
    class Meta:
        model = DownloaderJob
        fields = (
            "id",
            "downloader_task",
            "num_retries",
            "retried",
            "was_recreated",
            "worker_id",
            "worker_version",
            "nomad_job_id",
            "failure_reason",
            "success",
            "original_files",
            "start_time",
            "end_time",
            "created_at",
            "last_modified",
        )


class ProcessorJobSerializer(serializers.ModelSerializer):
    class Meta:
        model = ProcessorJob
        fields = (
            "id",
            "pipeline_applied",
            "num_retries",
            "retried",
            "worker_id",
            "ram_amount",
            "volume_index",
            "worker_version",
            "failure_reason",
            "nomad_job_id",
            "success",
            "original_files",
            "datasets",
            "start_time",
            "end_time",
            "created_at",
            "last_modified",
        )


##
# Datasets
##


def validate_dataset(data):
    """ Basic dataset validation. Currently only checks formatting, not values. """
    if data["data"] is not None:
        if type(data["data"]) != dict:
            raise serializers.ValidationError("`data` must be a dict of lists.")

        for key, value in data["data"].items():
            if type(value) != list:
                raise serializers.ValidationError(
                    "`data` must be a dict of lists. Problem with `" + str(key) + "`"
                )

            try:
                if len(value) != len(set(value)):
                    raise serializers.ValidationError("Duplicate values detected in " + str(value))
            except Exception as e:
                raise serializers.ValidationError("Received bad dataset data: " + str(e))

    else:
        raise serializers.ValidationError("`data` must be a dict of lists.")


class CreateDatasetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dataset
        fields = ("id", "data", "email_address", "email_ccdl_ok")

    def validate(self, data):
        """
        Ensure this is something we want in our dataset.
        """
        try:
            validate_dataset(data)
        except Exception:
            raise
        return data


class DatasetDetailsExperimentSerializer(serializers.ModelSerializer):
    """ This serializer contains all of the information about an experiment needed for the download
    page
    """

    organisms = serializers.ReadOnlyField(source="organism_names")
    sample_metadata = serializers.ReadOnlyField(source="sample_metadata_fields")

    class Meta:
        model = Experiment
        fields = ("title", "accession_code", "organisms", "sample_metadata", "technology")


class DatasetSerializer(serializers.ModelSerializer):
    start = serializers.NullBooleanField(required=False)
    experiments = DatasetDetailsExperimentSerializer(
        source="get_experiments", many=True, read_only=True
    )
    organism_samples = serializers.SerializerMethodField(read_only=True)
    worker_version = serializers.SerializerMethodField(read_only=True)

    def __init__(self, *args, **kwargs):
        super(DatasetSerializer, self).__init__(*args, **kwargs)

        if "context" in kwargs:
            if "request" in kwargs["context"]:
                # only include the fields `experiments` and `organism_samples` when the param `?details=true`
                # is provided. This is used on the frontend to render the downloads page
                # thanks to https://django.cowhite.com/blog/dynamically-includeexclude-fields-to-django-rest-framwork-serializers-based-on-user-requests/
                if "details" not in kwargs["context"]["request"].query_params:
                    self.fields.pop("experiments")
                    self.fields.pop("organism_samples")
                    self.fields.pop("worker_version")

            # only include the field `download_url` if a valid token is specified
            # the token lookup happens in the view.
            if "token" not in kwargs["context"]:
                self.fields.pop("download_url")

    class Meta:
        model = Dataset
        fields = (
            "id",
            "data",
            "aggregate_by",
            "scale_by",
            "is_processing",
            "is_processed",
            "is_available",
            "has_email",
            "expires_on",
            "s3_bucket",
            "s3_key",
            "success",
            "failure_reason",
            "created_at",
            "last_modified",
            "start",
            "size_in_bytes",
            "sha1",
            "experiments",
            "organism_samples",
            "download_url",
            "quantile_normalize",
            "quant_sf_only",
            "svd_algorithm",
            "worker_version",
        )
        extra_kwargs = {
            "data": {"required": True,},
            "id": {"read_only": True,},
            "is_processing": {"read_only": True,},
            "is_processed": {"read_only": True,},
            "is_available": {"read_only": True,},
            "expires_on": {"read_only": True,},
            "s3_bucket": {"read_only": True,},
            "s3_key": {"read_only": True,},
            "success": {"read_only": True,},
            "failure_reason": {"read_only": True,},
            "created_at": {"read_only": True,},
            "last_modified": {"read_only": True,},
            "size_in_bytes": {"read_only": True,},
            "sha1": {"read_only": True,},
            "download_url": {"read_only": True,},
            "worker_version": {
                "read_only": True,
                "help_text": "Returns the latest version of refine.bio that was used to build this dataset.",
            },
        }

    def validate(self, data):
        """
        Ensure this is something we want in our dataset.
        """
        try:
            validate_dataset(data)
        except Exception:
            raise
        return data

    def get_organism_samples(self, obj):
        """
        Groups the sample accession codes inside a dataset by their organisms, eg:
        { HOMO_SAPIENS: [S1, S2], DANIO: [S3] }
        Useful to avoid sending sample information on the downloads page
        """
        samples = (
            obj.get_samples()
            .prefetch_related("organism")
            .values("organism__name", "accession_code")
            .order_by("organism__name", "accession_code")
        )

        result = defaultdict(list)
        for sample in samples:
            result[sample["organism__name"]].append(sample["accession_code"])

        return result

    def get_worker_version(self, obj):
        processor_jobs = obj.processor_jobs.order_by("-created_at").values_list(
            "worker_version", flat=True
        )
        if processor_jobs:
            return processor_jobs[0]
        else:
            return None


class APITokenSerializer(serializers.ModelSerializer):
    class Meta:
        model = APIToken
        fields = ("id", "is_activated", "terms_and_conditions")
        extra_kwargs = {
            "id": {"read_only": True},
            "is_activated": {"read_only": False},
            "terms_and_conditions": {"read_only": True},
        }


class CompendiumResultSerializer(serializers.ModelSerializer):
    primary_organism = serializers.StringRelatedField(read_only=True)
    organisms = serializers.StringRelatedField(many=True, read_only=True)
    computed_file = ComputedFileSerializer(source="get_computed_file", read_only=True)

    class Meta:
        model = CompendiumResult
        fields = (
            "id",
            "primary_organism",
            "organisms",
            "svd_algorithm",
            "quant_sf_only",
            "compendium_version",
            "computed_file",
        )
        read_only_fields = fields


class CompendiumResultWithUrlSerializer(serializers.ModelSerializer):
    primary_organism = serializers.StringRelatedField(read_only=True)
    organisms = serializers.StringRelatedField(many=True, read_only=True)
    computed_file = ComputedFileWithUrlSerializer(source="get_computed_file", read_only=True)

    class Meta:
        model = CompendiumResult
        fields = (
            "id",
            "primary_organism",
            "organisms",
            "svd_algorithm",
            "quant_sf_only",
            "compendium_version",
            "computed_file",
        )
        read_only_fields = fields


##
# ElasticSearch Document Serializers
##


class ExperimentDocumentSerializer(DocumentSerializer):
    """Serializer for the Experiment document."""

    class Meta(object):
        """Meta options."""

        document = ExperimentDocument
        fields = (
            "id",
            "title",
            "publication_title",
            "description",
            "technology",
            "accession_code",
            "alternate_accession_code",
            "submitter_institution",
            "has_publication",
            "publication_doi",
            "publication_authors",
            "sample_metadata_fields",
            "platform_names",
            "platform_accession_codes",
            "organism_names",
            "pubmed_id",
            "num_total_samples",
            "num_processed_samples",
            "num_downloadable_samples",
            "source_first_published",
        )
        read_only_fields = fields
