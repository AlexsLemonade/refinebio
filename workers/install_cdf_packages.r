library("AnnotationDbi")

source("https://bioconductor.org/biocLite.R")
biocLite("org.Bt.eg.db")
biocLite("org.Cf.eg.db")
biocLite("org.Gg.eg.db")

# The next few lines install the packages necessary for handling entrez ids.

packages <- row.names(available.packages(contriburl="http://brainarray.mbni.med.umich.edu/bioc/src/contrib/"))
#this fetches the files from the website

entrez_packages <- packages[grep("entrezg", packages)]
#this narrows down from all the files gotten just having entrezg in the name

install.packages(entrez_packages, repos="http://brainarray.mbni.med.umich.edu/bioc/", type="source")
#this installs the subset of the list from the previous step on your computer.
