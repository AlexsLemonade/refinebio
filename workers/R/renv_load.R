options(warn = 2)
options(repos = structure(c(
    Bioconductor = "https://bioconductor.org/packages/3.6/bioc",
    BioconductorAnnotation = "https://bioconductor.org/packages/3.6/data/annotation",
    BioconductorExperiment = "https://bioconductor.org/packages/3.6/data/experiment",
    CRAN = "https://cloud.r-project.org"
)))
options(Ncpus = parallel::detectCores())
options(renv.r.version = "3.4.4")
options(renv.settings.use.cache = FALSE)

install.packages("BiocInstaller")
install.packages("https://cloud.r-project.org/src/contrib/renv_0.16.0.tar.gz")

renv::consent(provided = TRUE)
renv::restore(prompt = FALSE, rebuild = TRUE)
