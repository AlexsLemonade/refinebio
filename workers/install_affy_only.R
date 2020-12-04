# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(repos=structure(c(CRAN="https://cloud.r-project.org")))
options(Ncpus=parallel::detectCores())

# Use devtools::install_version() to install packages in cran.
devtools::install_version('dplyr', version='1.0.0')
devtools::install_version('tidyr', version='1.1.0')
devtools::install_version('ff', version='2.2-14')

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
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/oligo_1.52.1.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/GEOquery_2.56.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/SCAN.UPC_2.30.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/affy_1.66.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/affyio_1.58.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/AnnotationDbi_1.52.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/zlibbioc_1.34.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/preprocessCore_1.50.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/genefilter_1.70.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/sva_3.36.0.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/tximport_1.16.1.tar.gz',
  'https://bioconductor.org/packages/3.11/bioc/src/contrib/limma_3.44.3.tar.gz'
)
install_with_url(bioc_pkgs)
