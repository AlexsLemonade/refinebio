# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(Ncpus=parallel::detectCores())
options(repos=structure(c(CRAN="https://cloud.r-project.org")))

devtools::install_url('https://bioconductor.org/packages/release/bioc/src/contrib/tximport_1.16.1.tar.gz')
