# bin/lib — rbio command modules

This package backs the `bin/rbio` entry script. Each user-facing
namespace lives in its own module; shared infrastructure lives in
underscored modules.

## Layout

```
_runtime.py     Globals, stderr, run, REPO_ROOT
_docker.py      bake_target, container checks, require_*
<namespace>.py  one per command namespace (build, db, es, test, ...)
```

Underscored modules are shared infrastructure imported by the
namespace modules. They don't appear in `rbio --help`.

## Registry

Each namespace module exports a single `COMMANDS` list of
`(name, fn, one-line help)` tuples:

```python
COMMANDS = [
    ("db:init",  cmd_db_init,  "bootstrap the data_refinery database + role"),
    ("db:psql",  cmd_db_psql,  "open a psql shell against data_refinery"),
    ("db:reset", cmd_db_reset, "drop + recreate data_refinery, then re-apply migrations"),
]
```

`bin/rbio` concatenates them into one `COMMAND_REGISTRY` and derives
both the dispatch dict and the top-level `--help` block from it.

## Adding a namespace

1. Create `bin/lib/<name>.py` with `cmd_*` functions and a `COMMANDS`
   list (see above for the tuple shape).
2. Add the module to `NAMESPACES` in `bin/rbio`.
