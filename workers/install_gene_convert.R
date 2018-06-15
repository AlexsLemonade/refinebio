options(repos=structure(c(CRAN="https://cloud.r-project.org")))

# Install dev packages
install.packages("devtools")

devtools::install_version('data.table', version='1.11.0')
devtools::install_version('optparse', version='1.4.4')

devtools::install_version('tidyverse', version='1.2.1')
devtools::install_version('rlang', version='0.2.1')
