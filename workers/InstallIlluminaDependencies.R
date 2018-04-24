install.packages("optparse")

install.packages("doParallel", repos="https://rweb.crmda.ku.edu/cran/")
install.packages("data.table", repos="https://rweb.crmda.ku.edu/cran/")

source("https://bioconductor.org/biocLite.R")
biocLite("oligo")
biocLite("limma")
biocLite("illuminaHumanv2.db")
biocLite("illuminaHumanv4.db")
