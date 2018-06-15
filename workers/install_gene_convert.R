options(repos=structure(c(CRAN="https://cloud.r-project.org")))

# Install dev packages
install.packages("devtools")
install.packages("dplyr")

devtools::install_version('data.table', version='1.11.0')
devtools::install_version('optparse', version='1.4.4')

install.packages("tidyverse")
install.packages("rlang")
