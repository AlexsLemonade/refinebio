# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(repos=structure(c(CRAN="https://cran.microsoft.com/snapshot/2019-07-03")))
options(Ncpus=parallel::detectCores())

# Bioconductor packages, installed by devtools::install_url()

# devtools::install_url() requires BiocInstaller
install.packages('https://bioconductor.org/packages/3.6/bioc/src/contrib/BiocInstaller_1.28.0.tar.gz')

# Helper function that installs a list of packages based on input URL
install_with_url <- function(main_url, packages) {
  lapply(packages,
         function(pkg) devtools::install_url(paste0(main_url, pkg)))
}

bioc_url <- 'https://bioconductor.org/packages/3.6/bioc/src/contrib/'
bioc_pkgs <- c(
  'oligo_1.42.0.tar.gz',
  'Biobase_2.38.0.tar.gz',
  'GEOquery_2.46.15.tar.gz',
  'SCAN.UPC_2.20.0.tar.gz',
  'affy_1.56.0.tar.gz',
  'affyio_1.48.0.tar.gz',
  'AnnotationDbi_1.40.0.tar.gz',
  'zlibbioc_1.24.0.tar.gz',
  'preprocessCore_1.40.0.tar.gz',
  'genefilter_1.60.0.tar.gz',
  'sva_3.26.0.tar.gz',
  'tximport_1.6.0.tar.gz',
  'limma_3.34.9.tar.gz'
)
install_with_url(bioc_url, bioc_pkgs)
