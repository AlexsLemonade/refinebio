# Data Refinery Workers

This is the project root for the Data Refinery Workers. This project is
composed of a number of Nomad jobs which can be used to download and
process data from a variety of sources.

## Developing

When developing a new task you will probably need to run the task repeatedly.
This can be done easily by running the workers with `./run_workers` and then
modifying the
`data_refinery_workers/downloaders/management/commands/queue_task.py` file
to run the task you're developing. Once you've done that you can queue the task
with `./run_tester.py`

The worker container is run with a name of `worker1` so that it's output can
easily be inspected with `docker logs worker1`. However this means that you
cannot run `./run_worker.sh` twice in a row without deleting the old container.
This can be done easily with
```
docker stop worker1 && docker container prune -f
```

A development workflow might look like:
```
./run_worker.sh
./run_tester.sh
docker logs worker1
# Review the output and make changes
docker stop worker1 && docker container prune -f
```
