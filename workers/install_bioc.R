# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)
options(repos=structure(c(CRAN="https://cloud.r-project.org")))
options(Ncpus=parallel::detectCores())

# Use devtools::install_version() to install packages in cran.
devtools::install_version('dplyr', version='1.0.0')
devtools::install_version('tidyr', version='1.1.0')

# devtools::install_url() requires BiocInstaller
# install.packages('https://bioconductor.org/packages/3.6/bioc/src/contrib/BiocInstaller_1.28.0.tar.gz')
