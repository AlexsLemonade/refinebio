# Data Refinery Foreman

This is the project root for the Data Refinery Foreman. This project is
responsible for surveying a list of sources to discover new data,
scheduling those surveying tasks, queuing jobs for downloading and processing,
and managing those jobs.

## Development

### Running the surveyor

To run a surveyor job in a Docker container:

```sh
rbio compose:manage foreman survey_all --accession <ACCESSION>
```

`compose:manage` runs any Django management command via
`docker compose run --rm foreman python3 manage.py <cmd> ...`.

You can also start an interactive python interpreter by running
`rbio compose:manage foreman shell`. The interpreter runs inside the
foreman container with the environment configured appropriately.

### Testing

To run all foreman tests: `rbio test:foreman`. Pass `-t <tag>` to filter
by Django test tag, or any other `manage.py test` argument.
