options(warn = 2)
options(repos = structure(c(
    Bioconductor = "https://bioconductor.org/packages/3.11/bioc/",
    BioconductorAnnotation = "https://bioconductor.org/packages/3.11/data/annotation/",
    BioconductorExperiment = "https://bioconductor.org/packages/3.11/data/experiment",
    CRAN = "https://cloud.r-project.org/"
)))
options(Ncpus = parallel::detectCores())
options(renv.settings.use.cache = FALSE)

install.packages("renv")

renv::init(bare = TRUE)
