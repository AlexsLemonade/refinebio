"""Shared helpers used by every test:* command."""

import shlex


def coverage_command(extra_args):
    """Bash one-liner: run manage.py under coverage, write XML, print report, preserve test exit code."""
    args_str = " ".join(shlex.quote(a) for a in extra_args)
    # chmod opens up the bind-mount output so the host (CI runner, dev user)
    # can read coverage.xml regardless of which UID the container runs as.
    return (
        f'coverage run --source="." manage.py test --settings=tests.settings --no-input {args_str}; '
        "exit_code=$?; "
        "coverage xml -o data_store/coverage.xml; "
        "chmod -R a+rw data_store; "
        "coverage report -m; "
        "exit $exit_code"
    )
