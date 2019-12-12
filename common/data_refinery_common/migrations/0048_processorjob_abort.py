from django.db import migrations, models


def batch_update(queryset, batch=1000, **changes):
    """ Update per batch """
    is_done = False
    current = 0

    while not is_done:
        current = current + 1
        start = (current - 1) * batch
        end = start + batch
        pks = queryset.values_list('pk', flat=True)[start:end][::1]
        if pks:
            queryset.model.objects.filter(pk__in=pks).update(**changes)

        is_done = len(pks) < batch


def set_default_abort_value(apps, schema_editor):
    ProcessorJob = apps.get_model('data_refinery_common', 'ProcessorJob')
    queryset = ProcessorJob.objects.all()
    batch_update(queryset, abort=False)


class Migration(migrations.Migration):

    dependencies = [
        ('data_refinery_common', '0047_clean_compendia_files'),
    ]

    # we need to add a new field `abort` with default value `False`
    # The processor jobs table is big, so it the migration might time out
    # that's why we add the default value in batches
    operations = [
        # 1. remove indexes
        migrations.RemoveIndex(
            model_name='processorjob',
            name='processor_jobs_created_at',
        ),
        # 2. add fields with None as default value
        migrations.AddField(
            model_name='processorjob',
            name='abort',
            field=models.BooleanField(null=True),
        ),
        # 3. set default value in batches
        migrations.RunPython(set_default_abort_value),
        # 4. add default value for field
        migrations.AlterField(
            model_name='processorjob',
            name='abort',
            field=models.BooleanField(default=False, null=False),
        ),
        # 5. re-add indexes
        migrations.AddIndex(
            model_name='processorjob',
            index=models.Index(fields=['created_at'], name='processor_jobs_created_at'),
        ),
    ]
