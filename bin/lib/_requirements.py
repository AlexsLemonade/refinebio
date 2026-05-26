"""@requires decorator and epilog helper for rbio commands.

Commands declare env vars + external tools they need via @requires. Future
debug:env and debug:deps commands will aggregate this registry to render
per-command and cross-command views. Per-command help text pulls the same
metadata via epilog_from().
"""


def requires(*, env=None, env_optional=None, tools=None):
    """Declare a command's env var + external tool requirements.

    env             — required env vars (debug:env flags them red when unset)
    env_optional    — {name: default} env vars (default shown in help)
    tools           — external commands the implementation calls (docker, aws, ...)
    """

    def decorator(fn):
        fn._requirements = {
            "env": list(env or []),
            "env_optional": dict(env_optional or {}),
            "tools": list(tools or []),
        }
        return fn

    return decorator


def epilog_from(fn, extra=""):
    """Build a help-text epilog block from a command's @requires metadata."""
    req = getattr(fn, "_requirements", {})
    lines = []
    if req.get("env"):
        lines.append("required env:")
        lines.extend(f"  {v}" for v in req["env"])
    if req.get("env_optional"):
        lines.append("optional env (defaults shown):")
        lines.extend(f"  {v}={d}" for v, d in req["env_optional"].items())
    if req.get("tools"):
        lines.append(f"requires: {', '.join(req['tools'])}")
    if extra:
        if lines:
            lines.append("")
        lines.append(extra)
    return "\n".join(lines)
