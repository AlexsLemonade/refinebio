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

install_package_version("jsonlite", "1.5")
install_package_version("mime", "0.6")
install_package_version("curl", "3.2")
install_package_version("openssl", "1.0.2")
install_package_version("R6", "2.3.0")
install_package_version("httr", "1.3.1")
install_package_version("digest", "0.6.18")
install_package_version("memoise", "1.1.0")
install_package_version("whisker", "0.3-2")
install_package_version("rstudioapi", "0.8")
install_package_version("git2r", "0.23.0")
install_package_version("withr", "2.1.2")
install_package_version("devtools", "1.13.6")
