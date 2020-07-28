# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(repos=structure(c(CRAN="https://cran.revolutionanalytics.com")))
options(Ncpus=parallel::detectCores())

# Helper function that installs a list of packages based on input URL
install_with_url <- function(main_url, packages) {
  lapply(packages,
         function(pkg) devtools::install_url(paste0(main_url, pkg)))
}

bioc_url <- 'https://bioconductor.org/packages/release/bioc/src/contrib/'
bioc_pkgs <- c(
  'oligo_1.52.0.tar.gz',
  'GEOquery_2.56.0.tar.gz',
  'SCAN.UPC_2.30.0.tar.gz',
  'affy_1.66.0.tar.gz',
  'affyio_1.58.0.tar.gz',
  'AnnotationDbi_1.50.3.tar.gz',
  'zlibbioc_1.34.0.tar.gz',
  'preprocessCore_1.50.0.tar.gz',
  'genefilter_1.70.0.tar.gz',
  'sva_3.36.0.tar.gz',
  'tximport_1.16.1.tar.gz',
  'limma_3.44.3.tar.gz'
)
install_with_url(bioc_url, bioc_pkgs)
