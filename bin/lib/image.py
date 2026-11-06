"""image:* — content-addressed hashes for docker image tags."""

import argparse
import hashlib
import json
import re
import subprocess

from lib._requirements import requires
from lib._runtime import REPO_ROOT, stderr

# ---------- shared primitives ----------


def _sha256_file(path):
    """sha256(path contents), hex. Matches `sha256sum`'s hash column."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _sha256sum_line(path, display_path):
    """One line matching `sha256sum`'s output: '<64hex>  <path>\\n'."""
    return f"{_sha256_file(path)}  {display_path}\n"


def _parse_upstream_image(dockerfile_rel):
    """Read `FROM <image>@${...}` from a Dockerfile; return the <image> ref."""
    text = (REPO_ROOT / dockerfile_rel).read_text()
    m = re.search(r"^FROM\s+(\S+?)@\$\{", text, re.MULTILINE)
    if not m:
        raise RuntimeError(f"{dockerfile_rel}'s FROM line doesn't match 'FROM <image>@${{...}}'")
    return m.group(1)


def _resolve_upstream_digest(image_ref):
    """Return the linux/amd64 digest of `image_ref` as 'sha256:...', or None."""
    try:
        output = subprocess.check_output(
            ["docker", "manifest", "inspect", image_ref],
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        stderr(
            f"docker manifest inspect {image_ref} failed:\n"
            f"  {e.stderr.decode().strip() if e.stderr else e}"
        )
        return None
    data = json.loads(output)
    manifests = data.get("manifests")
    if manifests:
        for m in manifests:
            plat = m.get("platform", {})
            if plat.get("architecture") == "amd64" and plat.get("os") == "linux":
                return m["digest"]
        stderr(f"no linux/amd64 variant in {image_ref} manifest list")
        return None
    stderr(f"{image_ref} returned a non-list manifest; cannot pin")
    return None


# ---------- base: R 3.4.4 on pinned ubuntu ----------

_BASE_DOCKERFILE = "common/dockerfiles/Dockerfile.base"
_BASE_FILES = (_BASE_DOCKERFILE, "setup.cfg")
_BASE_DIRS = ("config",)
_BASE_PREFIX = "r3.4.4-"


def _hash_base(upstream_ref, upstream_digest):
    chunks = []
    for rel in _BASE_FILES:
        chunks.append(_sha256sum_line(REPO_ROOT / rel, rel))
    # Mirror the shell recipe `find <dir> -exec sha256sum {} + | sort` — sort
    # is by whole line (hash-first), not path, so produce lines then sort bytes.
    for d in _BASE_DIRS:
        root = REPO_ROOT / d
        dir_lines = []
        for p in root.rglob("*"):
            if p.is_file():
                rel = p.relative_to(REPO_ROOT).as_posix()
                dir_lines.append(_sha256sum_line(p, rel))
        dir_lines.sort()
        chunks.extend(dir_lines)
    chunks.append(f"{upstream_digest}  {upstream_ref}\n")
    return hashlib.sha256("".join(chunks).encode()).hexdigest()[:12]


def resolve_base(ubuntu_digest=None):
    """Resolve dr_base's content hash + the ubuntu digest it was built against.
    Returns (hash, digest) or (None, None) if the digest can't be fetched."""
    upstream_ref = _parse_upstream_image(_BASE_DOCKERFILE)
    digest = ubuntu_digest or _resolve_upstream_digest(upstream_ref)
    if digest is None:
        return None, None
    return _hash_base(upstream_ref, digest), digest


# ---------- r_env: per-worker R dependency layers on top of base ----------

WORKER_R_ENV_NAMES = (
    "smasher",
    "illumina",
    "salmon",
    "no_op",
    "downloaders",
    "compendia",
    "affymetrix",
)


def _r_env_files(worker):
    """Inputs that invalidate one worker's dr_<worker>_r_env image."""
    if worker == "affymetrix":
        dockerfile = "workers/dockerfiles/Dockerfile.affymetrix_r_env"
        extra = ("workers/R/dependencies/affymetrix/install_ensg_pkgs.R",)
    else:
        dockerfile = "workers/dockerfiles/Dockerfile.r_env"
        extra = ()
    return (
        "common/R/renv_load.R",
        f"workers/R/dependencies/{worker}/renv.lock",
        dockerfile,
        *extra,
    )


def _hash_r_env(worker, base_hash):
    chunks = [f"{base_hash}  parent:base\n"]
    for rel in _r_env_files(worker):
        chunks.append(_sha256sum_line(REPO_ROOT / rel, rel))
    return hashlib.sha256("".join(chunks).encode()).hexdigest()[:12]


def resolve_worker_r_env(worker, base_hash=None):
    """Resolve one worker's r_env hash. Resolves dr_base too if base_hash
    isn't supplied. Returns None if the base can't be resolved."""
    if worker not in WORKER_R_ENV_NAMES:
        raise ValueError(f"unknown r_env worker: {worker}")
    if base_hash is None:
        base_hash, _ = resolve_base()
        if base_hash is None:
            return None
    return _hash_r_env(worker, base_hash)


# ---------- CLI ----------

_IMAGE_NAMES = ("base",) + tuple(f"{w}_r_env" for w in WORKER_R_ENV_NAMES)


@requires(tools=["docker"])
def cmd_image_hash(argv):
    image_list = f"  base            (prefix: {_BASE_PREFIX})\n" + "\n".join(
        f"  {w}_r_env" for w in WORKER_R_ENV_NAMES
    )
    p = argparse.ArgumentParser(
        prog="rbio image:hash",
        description="Compute the content-addressed tag suffix for a docker image.",
        epilog=(
            "hashable images:\n"
            f"{image_list}\n"
            "\n"
            "examples:\n"
            "  rbio image:hash base\n"
            "  rbio image:hash base --ubuntu-digest=sha256:abc...\n"
            "  rbio image:hash smasher_r_env\n"
            "  rbio image:hash smasher_r_env --base-hash=0418c52de898\n"
            "\n"
            "output:\n"
            "  prints the content-addressed tag suffix on stdout.\n"
            "  use as a docker tag, e.g. dr_base:$(rbio image:hash base)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "image", choices=_IMAGE_NAMES, help=f"image to hash (one of: {', '.join(_IMAGE_NAMES)})"
    )
    p.add_argument(
        "--ubuntu-digest",
        default=None,
        metavar="sha256:...",
        help="fold a specific ubuntu digest into the base hash (base only)",
    )
    p.add_argument(
        "--base-hash",
        default=None,
        metavar="HASH",
        help="override the parent BASE_HASH (r_env images only)",
    )
    args = p.parse_args(argv)

    try:
        if args.image == "base":
            if args.base_hash is not None:
                stderr("rbio image:hash: --base-hash is not valid for 'base'")
                return 2
            h, _ = resolve_base(ubuntu_digest=args.ubuntu_digest)
            if h is None:
                return 1
            print(_BASE_PREFIX + h)
            return 0

        if args.ubuntu_digest is not None:
            stderr(f"rbio image:hash: --ubuntu-digest only valid for 'base', not '{args.image}'")
            return 2
        worker = args.image.removesuffix("_r_env")
        h = resolve_worker_r_env(worker, base_hash=args.base_hash)
        if h is None:
            return 1
        print(h)
        return 0
    except FileNotFoundError as e:
        stderr(f"rbio image:hash: missing input file: {e.filename}")
        return 1


COMMANDS = [
    ("image:hash", cmd_image_hash, "compute the content-addressed tag suffix for a docker image"),
]
