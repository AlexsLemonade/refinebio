##
# Contains helper serializers needed by other serializers
# * These are the same as the old serializers from serializers.py used by other serializers except with "Relation" added in the name
# * Some of these relation serializers are dependent on each other
##

from rest_framework import serializers

from data_refinery_common.models import (
    ComputationalResult,
    ComputationalResultAnnotation,
    ComputedFile,
    DownloaderJob,
    Organism,
    OrganismIndex,
    Processor,
    ProcessorJob,
    Sample,
)

##
# Organisms
##


class OrganismRelationSerializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = (
            "name",
            "taxonomy_id",
        )


##
# Computed Files
##


class ComputedFileRelationSerializer(serializers.ModelSerializer):
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


class ComputedFileWithUrlRelationSerializer(serializers.ModelSerializer):
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


##
# Samples
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


class SampleRelationSerializer(serializers.ModelSerializer):
    organism = OrganismRelationSerializer(many=False)

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


##
# Jobs
##


class DownloaderJobRelationSerializer(serializers.ModelSerializer):
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


class ProcessorJobRelationSerializer(serializers.ModelSerializer):
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
# Transcriptome Index
##


class OrganismIndexRelationSerializer(serializers.ModelSerializer):

    organism_name = serializers.StringRelatedField(source="organism", read_only=True)
    download_url = serializers.SerializerMethodField()

    class Meta:
        model = OrganismIndex

        fields = (
            "id",
            "assembly_name",
            "organism_name",
            "database_name",
            "release_version",
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
# Processors
##


class ProcessorRelationSerializer(serializers.ModelSerializer):
    class Meta:
        model = Processor
        fields = ("id", "name", "version", "docker_image", "environment")


##
# Results
##


class ComputationalResultAnnotationRelationSerializer(serializers.ModelSerializer):
    class Meta:
        model = ComputationalResultAnnotation
        fields = ("id", "data", "is_ccdl", "created_at", "last_modified")


class ComputationalResultRelationSerializer(serializers.ModelSerializer):
    annotations = ComputationalResultAnnotationRelationSerializer(
        many=True, source="computationalresultannotation_set"
    )
    processor = ProcessorRelationSerializer(many=False)
    organism_index = OrganismIndexRelationSerializer(many=False)
    files = ComputedFileRelationSerializer(many=True, source="computedfile_set")

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


class ComputationalResultNoFilesRelationSerializer(serializers.ModelSerializer):
    annotations = ComputationalResultAnnotationRelationSerializer(
        many=True, source="computationalresultannotation_set"
    )
    processor = ProcessorRelationSerializer(many=False)
    organism_index = OrganismIndexRelationSerializer(many=False)

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
