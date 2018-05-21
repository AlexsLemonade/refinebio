options(repos=structure(c(CRAN="https://cloud.r-project.org")))

# Install dev packages
install.packages("devtools")

devtools::install_version('doParallel', version='1.0.11')
devtools::install_version('data.table', version='1.11.0')
devtools::install_version('optparse', version='1.4.4')

# devtools::install_url() requires biocLite.R
source('https://bioconductor.org/biocLite.R')

# Helper function that installs a list of packages based on input URL
install_with_url <- function(main_url, packages) {
  lapply(packages,
         function(pkg) devtools::install_url(paste0(main_url, pkg)))
}

bioc_url <- 'https://bioconductor.org/packages/3.6/bioc/src/contrib/'
bioc_pkgs <- c(
   'oligo_1.42.0.tar.gz',
   'limma_3.34.9.tar.gz'
)
install_with_url(bioc_url, bioc_pkgs)

release_url <- 'https://bioconductor.org/packages/release/data/annotation/src/contrib/'
illumina_pkgs <- c(
  'illuminaHumanv2.db_1.26.0.tar.gz',
  'illuminaHumanv4.db_1.26.0.tar.gz'
)
install_with_url(release_url, illumina_pkgs)
