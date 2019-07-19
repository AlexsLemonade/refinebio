options(warn=2)
options(repos=structure(c(CRAN="http://cran.us.r-project.org")))
options(Ncpus=parallel::detectCores())

devtools::install_version('doParallel', version='1.0.11')
devtools::install_version('data.table', version='1.11.0')
devtools::install_version('optparse', version='1.4.4')
devtools::install_version('lazyeval', version='0.2.1')
devtools::install_version('rlang', version='0.3.1')
devtools::install_version('vctrs', version='0.1.0')
devtools::install_version('blob', version='1.1.1')
devtools::install_version('pillar', version='1.3.1')
devtools::install_version('tibble', version='2.0.1')
devtools::install_version('dplyr', version='0.7.8')

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
   'limma_3.34.9.tar.gz',
   'AnnotationDbi_1.40.0.tar.gz'
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
library(AnnotationDbi)
library(lazyeval)
library(limma)
library(oligo)
library(doParallel)
