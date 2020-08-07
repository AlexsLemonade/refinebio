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

### Debugging GitHub Actions workflows

CircleCI had a nice button called "Rerun job with SSH" that made debugging
simple, but GitHub hasn't added that feature (yet). I was trying to debug our
end-to-end tests timing out, which meant I needed access to the nomad logs, so
I looked into alternative ways to SSH into the machine while it was running the
tests.

I tried an action called `action-tmate`, but I think `tmate` times out after a
certain amount of time, so I coundn't stay ssh-ed in long enough to actually
observe what was wrong.

What I settled on was adapted from [this
blogpost](https://dev.to/retyui/how-debugging-github-actions-with-ssh-273n):

1. Fork the refinebio repo (this lets you set your own secrets on it and run
   tests without bothering other people).

2. Create an [ngrok](https://ngrok.com/) account. ngrok is a service that lets
   you tunnel local connections to a public URL hosted by them, including SSH
   connections.

3. Add this line to `config.yaml` in the workflow you are trying to debug:

```yaml
- name: Start SSH via Ngrok
  run: curl -sL https://gist.githubusercontent.com/retyui/7115bb6acf151351a143ec8f96a7c561/raw/7099b9db76729dc5761da72aa8525f632d8875c9/debug-github-actions.sh | bash
  env:
    NGROK_TOKEN: ${{ secrets.NGROK_TOKEN }}
    USER_PASS: ${{ secrets.USER_PASS }}
```

4. Set the `USER_PASS` and `NGROK_TOKEN` secrets in your fork. `USER_PASS`
   should be the password you want to use when logging in over SSH.

5. Optionally add this line to the end of the workflow you are trying to debug.
   It keeps the runner alive for an hour after something fails:

```yaml
- name: Don't kill instace
  if: ${{ failure() }}
  run: sleep 1h # Prevent to killing instance after failure
```

Now, when you re-run the action the `Start SSH via Ngrok` step of the workflow
will print out a URL that you can SSH into using the password set in
`USER_PASS` and you can start debugging the failed workflow. Our repo is in the
`work/refinebio/refinebio` folder.
