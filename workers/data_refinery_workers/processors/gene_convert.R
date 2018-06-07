###
# gene_convert.R
# Alex's Lemonad Stand Foundation
# Childhood Cancer Data Lab
#
###

#######################
# The command interface!
#######################

library("optparse")
library(data.table)

option_list = list(
  make_option(c("-p", "--platform"), type="character", default="illuminaHumanv4", 
              help="Platform", metavar="character"),
  make_option(c("-i", "--inputFile"), type="character", default="", 
              help="inputFile", metavar="character"),
  make_option(c("-o", "--outputFile"), type="character", default="", 
              help="outputFile", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

platform <- opt$platform
filePath <- opt$inputFile
outFilePath <- opt$outputFile

# Get and use this DB
message("Loading...")
source("http://bioconductor.org/biocLite.R")

#biocLite(paste(platform, ".db", sep="")) # This is handled at the image level
library(paste(platform, ".db", sep=""), character.only=TRUE)
library(Biobase)
library("AnnotationDbi")

# Read the data file
message("Reading data file...")
suppressWarnings(data <- fread(filePath, stringsAsFactors=FALSE, sep="\t", header=TRUE, autostart=10, data.table=FALSE, check.names=FALSE, fill=TRUE, na.strings="", showProgress=FALSE))

# Build a PROBEID:ENSEMBL map
res <- select(get(paste0(platform, ".db", sep="")), keys=data$ID_REF, columns=c("ENSEMBL"), keytype="PROBEID")

# Replace the probe IDs with Ensembles
i <- 0
for(m in data[, 1]){
  data[i, 1] = res[i, 2]
  i <- i + 1;
}

# Remove all of the unmapped (NA) values (likely control probes?)
data <- data[complete.cases(data), ]

# Data here can have duplicate rows that need to be squished together (via mean/max/etc),
# but this should be done at smash-time so that we can provide options on the squish method.

# Save to output file
write.table(data, outFilePath, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
