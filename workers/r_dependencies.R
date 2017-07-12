install.packages('ff', repos='https://cloud.r-project.org/')
install.packages('XML', repos='https://cloud.r-project.org/')
install.packages('RCurl', repos='https://cloud.r-project.org/')

source('https://bioconductor.org/biocLite.R')
biocLite('oligo')
biocLite('Biobase')
biocLite('SCAN.UPC')
biocLite('affy')
biocLite('affyio')
biocLite('AnnotationDbi')
biocLite('org.Hs.eg.db')
biocLite('org.Mm.eg.db')
biocLite('org.Dm.eg.db')
biocLite('org.Ce.eg.db')
biocLite('org.Bt.eg.db')
biocLite('org.Cf.eg.db')
biocLite('org.Gg.eg.db')
biocLite('org.Rn.eg.db')
biocLite('org.Ss.eg.db')
biocLite('org.Dr.eg.db')


packages <- row.names(available.packages(contriburl="http://brainarray.mbni.med.umich.edu/bioc/src/contrib/"))

## This narrows down from all the files gotten just having entrezg in the name
entrez_packages <- packages[grep("entrezg", packages)]

install.packages(entrez_packages, repos="http://brainarray.mbni.med.umich.edu/bioc/", type="source")
