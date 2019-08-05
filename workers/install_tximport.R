# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(Ncpus=parallel::detectCores())
options(repos=structure(c(CRAN="https://cran.microsoft.com/snapshot/2019-07-03")))

devtools::install_version('optparse', version='1.4.4')
devtools::install_version('rjson', version='0.2.19')
devtools::install_version('readr', version='1.1.1')

# devtools::install_url() requires BiocInstaller
install.packages('https://bioconductor.org/packages/3.6/bioc/src/contrib/BiocInstaller_1.28.0.tar.gz')

devtools::install_url('https://bioconductor.org/packages/3.6/bioc/src/contrib/tximport_1.6.0.tar.gz')
