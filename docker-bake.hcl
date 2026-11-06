// docker-bake.hcl
//
// refinebio image build configuration for `docker buildx bake`.
//
// Usage (local, no Dockerhub creds):
//   docker buildx bake base
//   docker buildx bake default
//   docker buildx bake all
//
// Usage (CI / push):
//   DOCKERHUB_REPO=ccdlstaging SYSTEM_VERSION=$CI_TAG PUSH_CACHE=true \
//     docker buildx bake all --push

// ---------------------------------------------------------------------------
// Variables - all overridable via environment.
// ---------------------------------------------------------------------------

// Default matches the project's existing tag convention so local builds are
// findable by tag without a registry round-trip.
// TODO: the deploy box runs Ubuntu 18.04 with Buildx v0.10.5, which does not support
// the `description` field in variable blocks (requires Buildx v0.13+). Once the deploy
// box is upgraded to Ubuntu 22.04+, uncomment the description lines below.
variable "DOCKERHUB_REPO" {
  // description = "Registry namespace used for image tags (e.g. <repo>/dr_base:<version>)."
  default = "ccdlstaging"
}

variable "SYSTEM_VERSION" {
  // description = "Tag suffix applied to every built image and cache ref."
  default = "local"
}

// Content-addressed tag inputs. All default to empty; the caller must set
// them in the environment when building any content-hashed target,
// otherwise the FROM lines parse to invalid references.
variable "UBUNTU_DIGEST" {
  // description = "sha256 digest of ubuntu:22.04 to pin dr_base's FROM to."
  default = ""
}

variable "BASE_HASH" {
  // description = "Content hash of dr_base inputs; drives dr_base's :r3.4.4-<hash> tag."
  default = ""
}

// Per-worker content hashes. Each is independent; a change to one worker's
// renv.lock invalidates only that worker's r_env image.
variable "SMASHER_R_ENV_HASH"     { default = "" }
variable "ILLUMINA_R_ENV_HASH"    { default = "" }
variable "SALMON_R_ENV_HASH"      { default = "" }
variable "NO_OP_R_ENV_HASH"       { default = "" }
variable "DOWNLOADERS_R_ENV_HASH" { default = "" }
variable "COMPENDIA_R_ENV_HASH"   { default = "" }
variable "AFFYMETRIX_R_ENV_HASH"  { default = "" }

// Default is the public read-only cache that CI pushes to.
// Note: when DOCKERHUB_REPO differs, cache is also read from DOCKERHUB_REPO
// so writers benefit from their own previously-pushed cache.
variable "CACHE_FROM_REPO" {
  // description = "Registry to read cache layers from."
  default = "ccdlstaging"
}

// Default false so unauthed local builds don't fail trying to push cache.
variable "PUSH_CACHE" {
  // description = "If true, write cache layers to DOCKERHUB_REPO. Requires push access."
  default = false
}

// Default amd64 because the Salmon / SRA Toolkit binaries (Dockerfile.salmon,
// Dockerfile.transcriptome) are x86_64-only and some apt packages have no
// arm64 variant.
variable "PLATFORMS" {
  // description = "Comma-separated platforms to build for."
  default = "linux/amd64"
}

// ---------------------------------------------------------------------------
// Helper functions.
// ---------------------------------------------------------------------------

function "tags_for" {
  params = [name]
  result = [
    "${DOCKERHUB_REPO}/dr_${name}:${SYSTEM_VERSION}",
    "${DOCKERHUB_REPO}/dr_${name}:latest",
  ]
}

function "cache_from_for" {
  params = [name]
  result = concat(
    [
      "type=registry,ref=${CACHE_FROM_REPO}/dr_${name}_cache:latest",
      "type=registry,ref=${CACHE_FROM_REPO}/dr_${name}_cache:${SYSTEM_VERSION}",
    ],
    DOCKERHUB_REPO != CACHE_FROM_REPO ? [
      "type=registry,ref=${DOCKERHUB_REPO}/dr_${name}_cache:latest",
      "type=registry,ref=${DOCKERHUB_REPO}/dr_${name}_cache:${SYSTEM_VERSION}",
    ] : []
  )
}

function "cache_to_for" {
  params = [name]
  result = PUSH_CACHE ? [
    "type=registry,ref=${DOCKERHUB_REPO}/dr_${name}_cache:latest,mode=max",
    "type=registry,ref=${DOCKERHUB_REPO}/dr_${name}_cache:${SYSTEM_VERSION},mode=max",
  ] : []
}

// Wire a child target's `FROM` line to a sibling parent target.
// Eliminates the registry round-trip that
//   FROM ${DOCKERHUB_REPO}/dr_<parent>:${SYSTEM_VERSION}
// would otherwise need.
function "from_target" {
  params = [parent]
  result = {
    "${DOCKERHUB_REPO}/dr_${parent}:${SYSTEM_VERSION}" = "target:${parent}"
  }
}

// Same, for content-addressed parents whose tag is not SYSTEM_VERSION.
// Example: from_target_hashed("base", base_tag())
function "from_target_hashed" {
  params = [parent, tag]
  result = {
    "${DOCKERHUB_REPO}/dr_${parent}:${tag}" = "target:${parent}"
  }
}

// Single source of truth for the dr_base content-addressed tag. The r3.4.4
// prefix documents the R version baked into dr_base; r_env Dockerfiles use
// this tag (passed as BASE_TAG) rather than reconstructing it.
function "base_tag" {
  params = []
  result = "r3.4.4-${BASE_HASH}"
}

// ---------------------------------------------------------------------------
// Abstract base target - concrete targets inherit platform + build args.
// ---------------------------------------------------------------------------

target "_common" {
  context   = "."
  platforms = split(",", PLATFORMS)
  args = {
    DOCKERHUB_REPO = DOCKERHUB_REPO
    SYSTEM_VERSION = SYSTEM_VERSION
    UBUNTU_DIGEST  = UBUNTU_DIGEST
    BASE_TAG       = base_tag()
  }
}

// ---------------------------------------------------------------------------
// common/ images.
// ---------------------------------------------------------------------------

target "base" {
  inherits   = ["_common"]
  dockerfile = "common/dockerfiles/Dockerfile.base"
  // The content-addressed tag r3.4.4-<BASE_HASH> is what downstream content-
  // hashed targets (dr_<worker>_r_env) FROM. The legacy :$SYSTEM_VERSION and
  // :latest tags are kept so untouched downstream images still resolve.
  tags = concat(
    tags_for("base"),
    BASE_HASH != "" ? ["${DOCKERHUB_REPO}/dr_base:${base_tag()}"] : []
  )
  cache-from = cache_from_for("base")
  cache-to   = cache_to_for("base")
}

target "migrations" {
  inherits   = ["_common"]
  dockerfile = "common/dockerfiles/Dockerfile.migrations"
  contexts   = from_target("base")
  tags       = tags_for("migrations")
  cache-from = cache_from_for("migrations")
  cache-to   = cache_to_for("migrations")
}

target "common_tests" {
  inherits   = ["_common"]
  dockerfile = "common/dockerfiles/Dockerfile.common_tests"
  contexts   = from_target("base")
  tags       = tags_for("common_tests")
  cache-from = cache_from_for("common_tests")
  cache-to   = cache_to_for("common_tests")
}

// ---------------------------------------------------------------------------
// api/ images.
// ---------------------------------------------------------------------------

target "api_base" {
  inherits   = ["_common"]
  dockerfile = "api/dockerfiles/Dockerfile.api_base"
  tags       = tags_for("api_base")
  cache-from = cache_from_for("api_base")
  cache-to   = cache_to_for("api_base")
}

target "api" {
  inherits   = ["_common"]
  dockerfile = "api/dockerfiles/Dockerfile.api"
  contexts   = from_target("api_base")
  tags       = tags_for("api")
  cache-from = cache_from_for("api")
  cache-to   = cache_to_for("api")
}

target "api_local" {
  inherits   = ["_common"]
  dockerfile = "api/dockerfiles/Dockerfile.api_local"
  contexts   = from_target("api_base")
  tags       = tags_for("api_local")
  cache-from = cache_from_for("api_local")
  cache-to   = cache_to_for("api_local")
}

// ---------------------------------------------------------------------------
// foreman/.
// ---------------------------------------------------------------------------

target "foreman" {
  inherits   = ["_common"]
  dockerfile = "foreman/dockerfiles/Dockerfile.foreman"
  contexts   = from_target("base")
  tags       = tags_for("foreman")
  cache-from = cache_from_for("foreman")
  cache-to   = cache_to_for("foreman")
}

// ---------------------------------------------------------------------------
// workers/ images.
// ---------------------------------------------------------------------------

target "transcriptome" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.transcriptome"
  contexts   = from_target("base")
  tags       = tags_for("transcriptome")
  cache-from = cache_from_for("transcriptome")
  cache-to   = cache_to_for("transcriptome")
}

// Content-addressed intermediates: the slow renv restore per worker, tagged
// by content hash so they survive SYSTEM_VERSION bumps.

target "smasher_r_env" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.r_env"
  args       = { WORKER = "smasher" }
  contexts   = from_target_hashed("base", base_tag())
  tags       = SMASHER_R_ENV_HASH != "" ? ["${DOCKERHUB_REPO}/dr_smasher_r_env:${SMASHER_R_ENV_HASH}"] : []
  cache-from = cache_from_for("smasher_r_env")
  cache-to   = cache_to_for("smasher_r_env")
}

target "smasher" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.smasher"
  args       = { R_ENV_HASH = SMASHER_R_ENV_HASH }
  contexts   = from_target_hashed("smasher_r_env", SMASHER_R_ENV_HASH)
  tags       = tags_for("smasher")
  cache-from = cache_from_for("smasher")
  cache-to   = cache_to_for("smasher")
}

target "illumina_r_env" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.r_env"
  args       = { WORKER = "illumina" }
  contexts   = from_target_hashed("base", base_tag())
  tags       = ILLUMINA_R_ENV_HASH != "" ? ["${DOCKERHUB_REPO}/dr_illumina_r_env:${ILLUMINA_R_ENV_HASH}"] : []
  cache-from = cache_from_for("illumina_r_env")
  cache-to   = cache_to_for("illumina_r_env")
}

target "salmon_r_env" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.r_env"
  args       = { WORKER = "salmon" }
  contexts   = from_target_hashed("base", base_tag())
  tags       = SALMON_R_ENV_HASH != "" ? ["${DOCKERHUB_REPO}/dr_salmon_r_env:${SALMON_R_ENV_HASH}"] : []
  cache-from = cache_from_for("salmon_r_env")
  cache-to   = cache_to_for("salmon_r_env")
}

target "no_op_r_env" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.r_env"
  args       = { WORKER = "no_op" }
  contexts   = from_target_hashed("base", base_tag())
  tags       = NO_OP_R_ENV_HASH != "" ? ["${DOCKERHUB_REPO}/dr_no_op_r_env:${NO_OP_R_ENV_HASH}"] : []
  cache-from = cache_from_for("no_op_r_env")
  cache-to   = cache_to_for("no_op_r_env")
}

target "downloaders_r_env" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.r_env"
  args       = { WORKER = "downloaders" }
  contexts   = from_target_hashed("base", base_tag())
  tags       = DOWNLOADERS_R_ENV_HASH != "" ? ["${DOCKERHUB_REPO}/dr_downloaders_r_env:${DOWNLOADERS_R_ENV_HASH}"] : []
  cache-from = cache_from_for("downloaders_r_env")
  cache-to   = cache_to_for("downloaders_r_env")
}

target "compendia_r_env" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.r_env"
  args       = { WORKER = "compendia" }
  contexts   = from_target_hashed("base", base_tag())
  tags       = COMPENDIA_R_ENV_HASH != "" ? ["${DOCKERHUB_REPO}/dr_compendia_r_env:${COMPENDIA_R_ENV_HASH}"] : []
  cache-from = cache_from_for("compendia_r_env")
  cache-to   = cache_to_for("compendia_r_env")
}

target "affymetrix_r_env" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.affymetrix_r_env"
  contexts   = from_target_hashed("base", base_tag())
  tags       = AFFYMETRIX_R_ENV_HASH != "" ? ["${DOCKERHUB_REPO}/dr_affymetrix_r_env:${AFFYMETRIX_R_ENV_HASH}"] : []
  cache-from = cache_from_for("affymetrix_r_env")
  cache-to   = cache_to_for("affymetrix_r_env")
}

target "no_op" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.no_op"
  args       = { R_ENV_HASH = NO_OP_R_ENV_HASH }
  contexts   = from_target_hashed("no_op_r_env", NO_OP_R_ENV_HASH)
  tags       = tags_for("no_op")
  cache-from = cache_from_for("no_op")
  cache-to   = cache_to_for("no_op")
}

target "illumina" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.illumina"
  args       = { R_ENV_HASH = ILLUMINA_R_ENV_HASH }
  contexts   = from_target_hashed("illumina_r_env", ILLUMINA_R_ENV_HASH)
  tags       = tags_for("illumina")
  cache-from = cache_from_for("illumina")
  cache-to   = cache_to_for("illumina")
}

target "salmon" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.salmon"
  args       = { R_ENV_HASH = SALMON_R_ENV_HASH }
  contexts   = from_target_hashed("salmon_r_env", SALMON_R_ENV_HASH)
  tags       = tags_for("salmon")
  cache-from = cache_from_for("salmon")
  cache-to   = cache_to_for("salmon")
}

target "downloaders" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.downloaders"
  args       = { R_ENV_HASH = DOWNLOADERS_R_ENV_HASH }
  contexts   = from_target_hashed("downloaders_r_env", DOWNLOADERS_R_ENV_HASH)
  tags       = tags_for("downloaders")
  cache-from = cache_from_for("downloaders")
  cache-to   = cache_to_for("downloaders")
}

target "compendia" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.compendia"
  args       = { R_ENV_HASH = COMPENDIA_R_ENV_HASH }
  contexts   = from_target_hashed("compendia_r_env", COMPENDIA_R_ENV_HASH)
  tags       = tags_for("compendia")
  cache-from = cache_from_for("compendia")
  cache-to   = cache_to_for("compendia")
}

// Affymetrix is opt-in (heavy build, ~75 GB disk required).
target "affymetrix" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.affymetrix"
  args       = { R_ENV_HASH = AFFYMETRIX_R_ENV_HASH }
  contexts   = from_target_hashed("affymetrix_r_env", AFFYMETRIX_R_ENV_HASH)
  tags       = tags_for("affymetrix")
  cache-from = cache_from_for("affymetrix")
  cache-to   = cache_to_for("affymetrix")
}

// ---------------------------------------------------------------------------
// Groups.
// ---------------------------------------------------------------------------

group "workers" {
  targets = [
    "transcriptome",
    "smasher",
    "salmon",
    "no_op",
    "illumina",
    "downloaders",
    "compendia",
  ]
}

// `default` is every image except affymetrix (which is opt-in via the
// `affymetrix` target or the `all` group).
group "default" {
  targets = [
    "base",
    "migrations",
    "common_tests",
    "foreman",
    "api_base",
    "api",
    "api_local",
    "transcriptome",
    "smasher",
    "salmon",
    "no_op",
    "illumina",
    "downloaders",
    "compendia",
  ]
}

// `all` includes affymetrix and content-addressed intermediates.
group "all" {
  targets = [
    "base",
    "migrations",
    "common_tests",
    "foreman",
    "api_base",
    "api",
    "api_local",
    "transcriptome",
    "smasher_r_env",
    "smasher",
    "illumina_r_env",
    "illumina",
    "salmon_r_env",
    "salmon",
    "no_op_r_env",
    "no_op",
    "downloaders_r_env",
    "downloaders",
    "compendia_r_env",
    "compendia",
    "affymetrix_r_env",
    "affymetrix",
  ]
}

// `deploy` is the set CI builds and pushes during a tagged deploy.
// Excludes test-only images (migrations, common_tests, api_local).
group "deploy" {
  targets = [
    "base",
    "api_base",
    "api",
    "foreman",
    "smasher",
    "compendia",
    "illumina",
    "affymetrix",
    "salmon",
    "transcriptome",
    "no_op",
    "downloaders",
  ]
}
