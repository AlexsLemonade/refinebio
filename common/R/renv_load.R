options(warn = 2)
# Bioconductor 3.6 (R 3.4.x era) is archived. The main bioconductor.org URLs
# return 302 redirects that R 3.4.4's url() connections do not follow, causing
# install.packages() to fail on the PACKAGES index check. Use the archive
# mirror directly instead.
BIOC_ARCHIVE <- "https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.6"
options(repos = structure(c(
    Bioconductor = paste0(BIOC_ARCHIVE, "/bioc"),
    BioconductorAnnotation = paste0(BIOC_ARCHIVE, "/data/annotation"),
    BioconductorExperiment = paste0(BIOC_ARCHIVE, "/data/experiment"),
    CRAN = "https://cloud.r-project.org"
)))
options(Ncpus = parallel::detectCores())
options(renv.r.version = "3.4.4")
options(renv.settings.use.cache = FALSE)

install.packages("BiocInstaller")
install.packages("https://cran.r-project.org/src/contrib/Archive/renv/renv_0.16.0.tar.gz")

renv::consent(provided = TRUE)
renv::restore(prompt = FALSE, rebuild = TRUE)
