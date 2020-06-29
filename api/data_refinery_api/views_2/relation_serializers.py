from rest_framework import serializers

from data_refinery_common.models import (
    Sample,
    ComputedFile,
    ComputationalResult,
    ComputationalResultAnnotation,
    Processor,
    Organism,
    OrganismIndex,
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
