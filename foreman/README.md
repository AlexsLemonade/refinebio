# Data Refinery Foreman

This is the project root for the Data Refinery Foreman. This project is
responsible for surveying a list of sources to discover new data,
scheduling those surveying tasks, queuing jobs for downloading and processing,
and managing those jobs.

## Development

### Running the surveyor

To run the surveyor in a Docker container just run `./run_surveyor.sh`.

You can also start an interactive python interpreter by running `../run_shell.sh`.
This interpreter will actually be running within the surveyor Docker container
with the environment configured appropriately.

### Testing

To run all tests for the project just run `./run_tests.sh`. This script will
accept any arguments that are valid for `python -m unittest`.
