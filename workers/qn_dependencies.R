options(warn=2)
options(repos=structure(c(CRAN="https://cran.revolutionanalytics.com")))
options(Ncpus=parallel::detectCores())

# devtools::install_url() requires BiocInstaller
install.packages('https://bioconductor.org/packages/3.6/bioc/src/contrib/BiocInstaller_1.28.0.tar.gz')

# Helper function that installs a list of packages based on input URL
install_with_url <- function(main_url, packages) {
  lapply(packages,
         function(pkg) devtools::install_url(paste0(main_url, pkg)))
}

bioc_url <- 'https://bioconductor.org/packages/3.7/bioc/src/contrib/'
bioc_pkgs <- c(
   'preprocessCore_1.42.0.tar.gz'
)
install_with_url(bioc_url, bioc_pkgs)

# Load these libraries because apparently just installing them isn't
# enough to verify that they have complementary versions.
library("optparse")
library(data.table)
library("preprocessCore")
