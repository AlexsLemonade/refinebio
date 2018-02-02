# Turn warnings into errors because biocLite throws warnings instead
# of error if it fails to install something.
options(warn=2)

install.packages('ff', repos='https://cloud.r-project.org/')
install.packages('XML', repos='https://cloud.r-project.org/')
install.packages('RCurl', repos='https://cloud.r-project.org/')
install.packages('RSQLite', repos='https://cloud.r-project.org/')
install.packages('tibble', repos='https://cloud.r-project.org/')
install.packages('xtable', repos='https://cloud.r-project.org/')
install.packages('pkgconfig', repos='https://cloud.r-project.org/')

source('https://bioconductor.org/biocLite.R')
bioconductor.packages <- c(
  'oligo',
  'Biobase',
  'SCAN.UPC',
  'affy',
  'affyio',
  'AnnotationDbi',
  'zlibbioc',
  'preprocessCore',
  'genefilter',
  'sva',
  'org.Hs.eg.db',
  'org.Mm.eg.db',
  'org.Dm.eg.db',
  'org.Ce.eg.db',
  'org.Bt.eg.db',
  'org.Cf.eg.db',
  'org.Gg.eg.db',
  'org.Rn.eg.db',
  'org.Ss.eg.db',
  'org.Dr.eg.db'
)

lapply(bioconductor.packages, biocLite)


packages <- row.names(available.packages(contriburl="http://brainarray.mbni.med.umich.edu/bioc/src/contrib/"))

## This narrows down from all the files gotten just having entrezg in the name
entrez.packages <- packages[grep("entrezg", packages)]

install.packages(entrez.packages, repos="http://brainarray.mbni.med.umich.edu/bioc/", type="source")

platform.designs <- c(
  "pd.hugene.1.0.st.v1",
)

lapply(platform.designs, biocLite)
