"""test:* — per-subproject test runners (api / common / foreman / workers) + test:all."""

from lib.test.all import COMMANDS as _ALL_COMMANDS
from lib.test.subprojects import COMMANDS as _SUBPROJECT_COMMANDS
from lib.test.workers import COMMANDS as _WORKER_COMMANDS

COMMANDS = [*_SUBPROJECT_COMMANDS, *_WORKER_COMMANDS, *_ALL_COMMANDS]
