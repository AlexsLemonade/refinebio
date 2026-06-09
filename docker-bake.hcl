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

// ---------------------------------------------------------------------------
// Abstract base target - concrete targets inherit platform + build args.
// ---------------------------------------------------------------------------

target "_common" {
  context   = "."
  platforms = split(",", PLATFORMS)
  args = {
    DOCKERHUB_REPO = DOCKERHUB_REPO
    SYSTEM_VERSION = SYSTEM_VERSION
  }
}

// ---------------------------------------------------------------------------
// common/ images.
// ---------------------------------------------------------------------------

target "base" {
  inherits   = ["_common"]
  dockerfile = "common/dockerfiles/Dockerfile.base"
  tags       = tags_for("base")
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

target "smasher" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.smasher"
  contexts   = from_target("base")
  tags       = tags_for("smasher")
  cache-from = cache_from_for("smasher")
  cache-to   = cache_to_for("smasher")
}

target "salmon" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.salmon"
  contexts   = from_target("base")
  tags       = tags_for("salmon")
  cache-from = cache_from_for("salmon")
  cache-to   = cache_to_for("salmon")
}

target "no_op" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.no_op"
  contexts   = from_target("base")
  tags       = tags_for("no_op")
  cache-from = cache_from_for("no_op")
  cache-to   = cache_to_for("no_op")
}

target "illumina" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.illumina"
  contexts   = from_target("base")
  tags       = tags_for("illumina")
  cache-from = cache_from_for("illumina")
  cache-to   = cache_to_for("illumina")
}

target "downloaders" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.downloaders"
  contexts   = from_target("base")
  tags       = tags_for("downloaders")
  cache-from = cache_from_for("downloaders")
  cache-to   = cache_to_for("downloaders")
}

target "compendia" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.compendia"
  contexts   = from_target("base")
  tags       = tags_for("compendia")
  cache-from = cache_from_for("compendia")
  cache-to   = cache_to_for("compendia")
}

// Affymetrix is opt-in (heavy build, ~75 GB disk required).
target "affymetrix" {
  inherits   = ["_common"]
  dockerfile = "workers/dockerfiles/Dockerfile.affymetrix"
  contexts   = from_target("base")
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

// `all` includes affymetrix.
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
    "smasher",
    "salmon",
    "no_op",
    "illumina",
    "downloaders",
    "compendia",
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
