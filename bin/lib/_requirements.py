"""@requires decorator and epilog helper for rbio commands.

Commands declare env vars + external tools they need via @requires. The
debug:env / debug:deps commands aggregate this registry to render
per-command and cross-command views. Per-command help text pulls the same
metadata via epilog_from().
"""


def requires(*, env=None, env_optional=None, env_alternatives=None, tools=None):
    """Declare a command's env var + external tool requirements.

    env               — required env vars (debug:env flags them red when unset)
    env_optional      — {name: default} env vars (default shown in help)
    env_alternatives  — list of (label, options) where each option is a list of
                        vars that together (AND) satisfy that option; outer list
                        of options is OR. Example for AWS auth:
                          env_alternatives=[
                              ("AWS auth", [
                                  ["AWS_PROFILE"],
                                  ["AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY"],
                              ]),
                          ]
    tools             — external commands the implementation calls (docker, aws, ...)
    """

    def decorator(fn):
        fn._requirements = {
            "env": list(env or []),
            "env_optional": dict(env_optional or {}),
            "env_alternatives": list(env_alternatives or []),
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
    if req.get("env_alternatives"):
        for label, options in req["env_alternatives"]:
            opts = " or ".join(" + ".join(opt) for opt in options)
            lines.append(f"{label} (need one): {opts}")
    if req.get("tools"):
        lines.append(f"requires: {', '.join(req['tools'])}")
    if extra:
        if lines:
            lines.append("")
        lines.append(extra)
    return "\n".join(lines)
