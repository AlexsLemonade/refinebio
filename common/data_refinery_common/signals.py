from django.db.models.signals import pre_delete
from django.dispatch import receiver

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_common.models import ComputedFile, OriginalFile

logger = get_and_configure_logger(__name__)


@receiver(pre_delete, sender=OriginalFile)
@receiver(pre_delete, sender=ComputedFile)
def remove_file_from_s3(sender, instance, **kwargs):
    """ When the local model is about to be deleted, try to delete the s3 file
    """
    logger.info("pre_delete for sender {} called".format(sender))
    instance.delete_s3_file()
