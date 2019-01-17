"""This command will go through ComputationalResult objects that are a
result of the tximport Processor and associate them with the
experiment that tximport was run on. It's purpose is to populate
missing links in our data model, however once it does so there are
additional PRs planned which will break an assumption that this
command makes, so it should be removed once it has served its purpose.

The assumption this command is relying on is:
  * tximport has only been run on full experiments.
"""

from django.core.management.base import BaseCommand
from data_refinery_common.models import *


def make_experiment_result_associations():
        """This function performs the function explained at the head of this file.

        It does so by following this general strategy:
        1. Get tximport results by querying based on the processor
        2. Then get one associated sample
        3. Then go through that sample's experiments until an experiment is found
        that has all of its samples associated with that result.
        4. Then make an association with that result."""

        # There are multiple "Processor" objects for tximport because
        # we create a new one for each version. However we don't care
        # which version was used, we just all tximport results.
        tximport_processors = Processor.objects.filter(name="Tximport").all()

        for tximport_processor in tximport_processors:
            results = ComputationalResult.objects.filter(processor=tximport_processor)

            for result in results:
                result_sample = result.samples.first()

                for experiment in result_sample.experiments.all():
                    experiment_samples = experiment.samples.all()

                    num_result_associations = 0
                    for experiment_sample in experiment_samples:
                        try:
                            SampleResultAssociation.objects.get(sample=experiment_sample, result=result)

                            # If we've made it here, then the association exists so count it!
                            num_result_associations += 1
                        except:
                            # If we've made it here, then the
                            # association doesn't exist so this isn't
                            # the experiment that the result is for.
                            break

                    if num_result_associations == len(experiment_samples):
                        # Every sample in the experiment is associated
                        # with this ComputationalResult, so we can
                        # safely say the experiment is associated with
                        # it and make that relationship explicit.
                        ExperimentResultAssociation.objects.get_or_create(
                            experiment=experiment,
                            result=result
                        )


class Command(BaseCommand):
    def handle(self, *args, **options):
        """This is just the entrypoint for this management command.

        All of its work is done in a separate function because that
        makes it much easier to test."""
        make_experiment_result_associations()
