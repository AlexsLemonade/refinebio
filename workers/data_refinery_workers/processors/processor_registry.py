from data_refinery_workers.processors.array_express \
    import affy_to_pcl

"""
This is a dictionary which maps valid values for Batch.pipeline_required
to the processor pipeline Celery task.
"""

processor_pipeline_registry = {
    "AFFY_TO_PCL": affy_to_pcl
}
