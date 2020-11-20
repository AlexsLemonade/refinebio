# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(Ncpus=parallel::detectCores())
options(repos=structure(c(CRAN="https://cloud.r-project.org")))

# Helper function that installs a list of packages based on input URL
install_with_url <- function(main_url, packages) {
  pkg_ids <- devtools::install_url(paste0(main_url, packages))
  if(any(is.na(pkg_ids))) {
    pkg_fails <- paste(packages[is.na(pkg_ids)], collapse = "; ")
    stop(paste("Failed to install package(s):", pkg_fails ))
  }
  return(pkg_ids)
}

devtools::install_version('dplyr', version='1.0.2')

bioc_url <- 'https://bioconductor.org/packages/release/bioc/src/contrib/'
bioc_pkgs <- c(
  'AnnotationDbi_1.52.0.tar.gz'
)
install_with_url(bioc_url, bioc_pkgs)

release_url <- 'https://bioconductor.org/packages/release/data/annotation/src/contrib/'
illumina_pkgs <- c(
  'illuminaHumanv1.db_1.26.0.tar.gz',
  'illuminaHumanv2.db_1.26.0.tar.gz',
  'illuminaHumanv3.db_1.26.0.tar.gz',
  'illuminaHumanv4.db_1.26.0.tar.gz',
  'illuminaMousev1.db_1.26.0.tar.gz',
  'illuminaMousev1p1.db_1.26.0.tar.gz',
  'illuminaMousev2.db_1.26.0.tar.gz',
  'illuminaRatv1.db_1.26.0.tar.gz'
)
install_with_url(release_url, illumina_pkgs)

# Load these libraries because apparently just installing them isn't
# enough to verify that they have complementary versions.
library("optparse")
library(data.table)
library("dplyr")
library("rlang")
library(lazyeval)
library(AnnotationDbi)
