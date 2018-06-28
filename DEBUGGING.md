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
