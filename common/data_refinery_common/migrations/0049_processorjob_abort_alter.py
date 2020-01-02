from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("data_refinery_common", "0048_processorjob_abort"),
    ]

    # Final steps to add the abort field to the processor jobs table
    operations = [
        migrations.AlterField(
            model_name="processorjob",
            name="abort",
            field=models.BooleanField(default=False, null=False),
        ),
    ]
