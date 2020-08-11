### Debugging!

This file will serve as a collection of commands, queries, or other tricks for debugging our systems in production.
The hope is that it will grow organically as we debug various issues, however at some point as it grows it may need to be reorganized.
For now though, make each sub header the class of issue the section will help to debug.

#### Too Many ProcessorJobs Were Queued Per OriginalFile

In the scenario where for some reason too many ProcessorJobs were queued per OriginalFile, you may want to leave one ProcessorJob per OriginalFile so they will either run, or preserve the record of the work that was done.
The following queries in concert will leave one ProcessorJob per OriginalFile, and delete the rest.
The ProcessorJob with the smallest id will be left for reach OriginalFile.

WARNING: Any time you run a DELETE query, you should first replace `DELETE` with `SELECT *` or `SELECT COUNT(*)` to make sure that you know what you will be deleting and that it makes sense.


```SQL
DELETE FROM processor_jobs WHERE id NOT IN
(SELECT pj_id FROM
        (SELECT MIN(processorjob_originalfile_associations.processor_job_id) pj_id, original_file_id
         FROM processorjob_originalfile_associations
         GROUP BY original_file_id) AS pjs);
```

```SQL
DELETE FROM processorjob_originalfile_associations WHERE processor_job_id NOT IN
(SELECT pj_id FROM
        (SELECT MIN(processorjob_originalfile_associations.processor_job_id) pj_id, original_file_id
         FROM processorjob_originalfile_associations
         GROUP BY original_file_id) AS pjs);
```

### Emptying the queue

If you ever need to empty the queue of processor and downloader jobs, you can run this in the foreman:

```python
from data_refinery_common.models.jobs.processor_job import ProcessorJob
from data_refinery_common.models.jobs.downloader_job import DownloaderJob
from django.utils import timezone
from datetime import datetime

JOB_CREATED_AT_CUTOFF = datetime(2019, 9, 19, tzinfo=timezone.utc)

ProcessorJob.failed_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF).exclude(
    pipeline_applied="JANITOR"
).update(no_retry=True)
ProcessorJob.lost_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF).exclude(
    pipeline_applied="JANITOR"
).update(no_retry=True)
ProcessorJob.hung_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF).exclude(
    pipeline_applied="JANITOR"
).update(no_retry=True)

DownloaderJob.failed_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF).update(no_retry=True)
DownloaderJob.lost_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF).update(no_retry=True)
DownloaderJob.hung_objects.filter(created_at__gt=JOB_CREATED_AT_CUTOFF).update(no_retry=True)
```

This sets all lost, hung, and failed processor jobs to not retry.

### Summarizing failure reasons

If you want to figure out why processor jobs are failing, you can use this
command to summarize the failure reason for recent jobs (you can tweak what
"recent" means, here it means the last 10 hours).

```python
from data_refinery_common.models.jobs.processor_job import ProcessorJob
from django.utils import timezone
from datetime import timedelta
from collections import Counter
t = timezone.now() - timedelta(hours=10)

Counter([j.failure_reason for j in ProcessorJob.objects.filter(created_at__gt=t, success=False)])
```

or for the 5 most common failure reasons:

```python
Counter([j.failure_reason for j in ProcessorJob.objects.filter(created_at__gt=t, success=False)]).most_common(5)
```

You can also do the same for downloader jobs:

```python
from data_refinery_common.models.jobs.downloader_job import DownloaderJob
from django.utils import timezone
from datetime import timedelta
from collections import Counter
t = timezone.now() - timedelta(hours=10)

Counter([j.failure_reason for j in DownloaderJob.objects.filter(created_at__gt=t, success=False)])
```
