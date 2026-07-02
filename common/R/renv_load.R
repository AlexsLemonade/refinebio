options(warn = 2)
# Bioconductor 3.6 (R 3.4.x era) is archived. bioconductor.org's /packages/3.6
# URLs 302-redirect to the S3-hosted archive bucket, and R 3.4.4's url()
# connections do not follow redirects, so install.packages() fails on the
# PACKAGES index check. Point at the S3 bucket directly (path-style HTTPS so the
# dotted bucket name doesn't break TLS; no redirect for url() to follow).
BIOC_ARCHIVE <- "https://s3.amazonaws.com/archive.bioconductor.org/packages/3.6"
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
