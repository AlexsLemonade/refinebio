options(warn=2)
options(repos=structure(c(CRAN="https://cloud.r-project.org")))
options(Ncpus=parallel::detectCores())

# Helper function that installs a list of packages using the input URLs
install_with_url <- function(urls) {
  pkg_ids <- devtools::install_url(urls)
  if(any(is.na(pkg_ids))) {
    pkg_fails <- paste(urls[is.na(pkg_ids)], collapse = "; ")
    stop(paste("Failed to install package(s):", pkg_fails ))
  }
  return(pkg_ids)
}

devtools::install_version('dplyr', version='1.0.2')

bioc_pkgs <- c(
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/oligo_1.52.1.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/AnnotationDbi_1.52.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/limma_3.44.3.tar.gz'
)
install_with_url(bioc_pkgs)

illumina_pkgs <- c(
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaHumanv1.db_1.26.0.tar.gz',
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaHumanv2.db_1.26.0.tar.gz',
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaHumanv3.db_1.26.0.tar.gz',
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaHumanv4.db_1.26.0.tar.gz',
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaMousev1.db_1.26.0.tar.gz',
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaMousev1p1.db_1.26.0.tar.gz',
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaMousev2.db_1.26.0.tar.gz',
  'https://bioconductor.org/packages/3.11/data/annotation/src/contrib/illuminaRatv1.db_1.26.0.tar.gz'
)
install_with_url(illumina_pkgs)

# Load these libraries because apparently just installing them isn't
# enough to verify that they have complementary versions.
library("optparse")
library(data.table)
library("dplyr")
library("rlang")
library(AnnotationDbi)
library(lazyeval)
library(limma)
library(oligo)
library(doParallel)
