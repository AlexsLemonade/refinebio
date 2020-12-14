# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(repos=structure(c(CRAN="https://cloud.r-project.org")))
options(Ncpus=parallel::detectCores())

# Bioconductor packages, installed by devtools::install_url()


# Helper function that installs a list of packages using the input URLs
install_with_url <- function(urls) {
  pkg_ids <- devtools::install_url(urls)
  if(any(is.na(pkg_ids))) {
    pkg_fails <- paste(urls[is.na(pkg_ids)], collapse = "; ")
    stop(paste("Failed to install package(s):", pkg_fails ))
  }
  return(pkg_ids)
}

bioc_pkgs <- c(
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/affyio_1.58.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/zlibbioc_1.34.0.tar.gz'
)
install_with_url(bioc_pkgs)
