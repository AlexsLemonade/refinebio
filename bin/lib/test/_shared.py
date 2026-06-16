"""Shared helpers used by every test:* command."""

import shlex


def coverage_command(extra_args):
    """Bash one-liner: run manage.py under coverage, write XML, print report, preserve test exit code."""
    args_str = " ".join(shlex.quote(a) for a in extra_args)
    # Open up the bind-mount output so the host (CI runner, dev user) can read
    # coverage.xml regardless of which UID the container runs as. Filter to
    # files owned by the current UID — pre-existing committed fixtures belong
    # to the host user and chmod requires ownership, not just write bits.
    return (
        f'coverage run --source="." manage.py test --settings=tests.settings --no-input {args_str}; '
        "exit_code=$?; "
        "coverage xml -o data_store/coverage.xml; "
        'find data_store -user "$(id -u)" -exec chmod a+rw {} +; '
        "coverage report -m; "
        "exit $exit_code"
    )
