# This script installs version 1.13.6 of the R package devtools.

# Because devtools is the tool which will install specific versions of
# other R packages, it cannot be used it to install itself. However in
# order to ensure there are no dependency conflicts, it also needs to
# pin each version of its dependencies.

# This script uses install.packages to install all of the dependencies
# for devtools, however install.packages unfortunately does not
# recursively install dependencies for the packages it
# installs. Therefore this script installs every dependency up the
# chain to devtools. It's not ideal, but it guarantees a stable
# install of devtools that won't break when a new version of devtools
# is released.

# Cranlock was used to find the versions of dependencies to install

# Treat warnings as errors, set CRAN mirror, and set parallelization:
options(warn=2)
options(repos=structure(c(CRAN="https://cran.revolutionanalytics.com/")))
options(Ncpus=parallel::detectCores())


install_package_version <- function(package_name, version) {
  # This function install a specific version of a package.

  # However, because the most current version of a package lives in a
  # different location than the older versions, we have to check where
  # it can be found.
  package_tarball <- paste0(package_name, "_", version, ".tar.gz")
  package_url <- paste0("https://cran.revolutionanalytics.com/src/contrib/", package_tarball)

  # Give CRAN a full minute to timeout since it's not always the most reliable.
  curl_result <- system(paste0("curl --head --connect-timeout 60 ", package_url), intern=TRUE)
  if (grepl("404", curl_result[1])) {
    package_url <- paste0("https://cran.revolutionanalytics.com/src/contrib/Archive/", package_name, "/", package_tarball)

    # Make sure the package actually exists in the archive!
    curl_result <- system(paste0("curl --head --connect-timeout 120 ", package_url), intern=TRUE)
    if (grepl("404", curl_result[1])) {
      stop(paste("Package", package_name, "version", version, "does not exist!"))
    }
  }

  install.packages(package_url)
}

# Generated using cranlock
install_package_version("jsonlite", "1.7.0")
install_package_version("withr", "2.2.0")
install_package_version("crayon", "1.3.4")
install_package_version("fansi", "0.4.1")
install_package_version("glue", "1.4.1")
install_package_version("assertthat", "0.2.1")
install_package_version("cli", "2.0.2")
install_package_version("sessioninfo", "1.1.1")
install_package_version("rlang", "0.4.7")
install_package_version("ellipsis", "0.3.1")
install_package_version("xml2", "1.3.2")
install_package_version("curl", "4.3")
install_package_version("rversions", "2.0.2")
install_package_version("remotes", "2.2.0")
install_package_version("R6", "2.4.1")
install_package_version("praise", "1.0.0")
install_package_version("magrittr", "1.5")
install_package_version("digest", "0.6.25")
install_package_version("evaluate", "0.14")
install_package_version("backports", "1.1.8")
install_package_version("rprojroot", "1.3-2")
install_package_version("desc", "1.2.0")
install_package_version("prettyunits", "1.1.1")
install_package_version("ps", "1.3.3")
install_package_version("processx", "3.4.3")
install_package_version("callr", "3.4.3")
install_package_version("pkgbuild", "1.1.0")
install_package_version("rstudioapi", "0.11")
install_package_version("pkgload", "1.1.0")
install_package_version("testthat", "2.3.2")
install_package_version("mime", "0.9")
install_package_version("sys", "3.4")
install_package_version("askpass", "1.1")
install_package_version("openssl", "1.4.2")
install_package_version("httr", "1.4.2")
install_package_version("fs", "1.4.2")
install_package_version("git2r", "0.27.1")
install_package_version("whisker", "0.4")
install_package_version("clipr", "0.7.0")
install_package_version("ini", "0.3.1")
install_package_version("gh", "1.1.0")
install_package_version("utf8", "1.1.4")
install_package_version("lifecycle", "0.2.0")
install_package_version("vctrs", "0.3.2")
install_package_version("pillar", "1.4.6")
install_package_version("pkgconfig", "2.0.3")
install_package_version("tibble", "3.0.3")
install_package_version("rematch2", "2.1.2")
install_package_version("purrr", "0.3.4")
install_package_version("yaml", "2.2.1")
install_package_version("usethis", "1.6.1")
install_package_version("lazyeval", "0.2.2")
install_package_version("rex", "1.2.0")
install_package_version("covr", "3.5.0")
install_package_version("memoise", "1.1.0")
install_package_version("base64enc", "0.1-3")
install_package_version("htmltools", "0.5.0")
install_package_version("crosstalk", "1.1.0.1")
install_package_version("htmlwidgets", "1.5.1")
install_package_version("Rcpp", "1.0.5")
install_package_version("BH", "1.72.0-3")
install_package_version("later", "1.1.0.1")
install_package_version("promises", "1.1.1")
install_package_version("DT", "0.14")
install_package_version("commonmark", "1.7")
install_package_version("stringi", "1.4.6")
install_package_version("stringr", "1.4.0")
install_package_version("xfun", "0.16")
install_package_version("markdown", "1.1")
install_package_version("highr", "0.8")
install_package_version("knitr", "1.29")
install_package_version("brew", "1.0-6")
install_package_version("roxygen2", "7.1.1")
install_package_version("xopen", "1.0.0")
install_package_version("rcmdcheck", "1.3.3")
install_package_version("devtools", "2.3.1")
