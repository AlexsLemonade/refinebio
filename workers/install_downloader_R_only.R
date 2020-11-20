# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(repos=structure(c(CRAN="https://cloud.r-project.org")))
options(Ncpus=parallel::detectCores())

# Bioconductor packages, installed by devtools::install_url()

# Helper function that installs a list of packages based on input URL
install_with_url <- function(main_url, packages) {
  pkg_ids <- devtools::install_url(paste0(main_url, packages))
  if(any(is.na(pkg_ids))) {
    pkg_fails <- paste(packages[is.na(pkg_ids)], collapse = "; ")
    stop(paste("Failed to install package(s):", pkg_fails ))
  }
  return(pkg_ids)
}

bioc_url <- 'https://bioconductor.org/packages/release/bioc/src/contrib/'
bioc_pkgs <- c(
  'affyio_1.60.0.tar.gz',
  'zlibbioc_1.36.0.tar.gz'
)
install_with_url(bioc_url, bioc_pkgs)
