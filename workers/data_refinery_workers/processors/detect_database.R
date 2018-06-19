###
# detect_database.R
# Alex's Lemonad Stand Foundation
# Childhood Cancer Data Lab
#
###

##
# get & read in non-normalized supplementary data
# IDs = identifiers used by submitter 
# 
# for each version of db package for organism (e.g., v1-v4 in human)
#   calculate the % of IDs in the probes for current package
#   calculate the % of probes in the IDs for current package
# 
# save the calculated overlaps as SampleAnnotation (will be number of versions x 2)
# highest_overlap = which version has the highest % of IDs in the probes
# 
# if highest_overlap > some high threshold
#   supply as platform/package to Illumina SCAN Rscript
# else
#   we should get the processed data and attempt to convert the IDs to Ensembl gene IDs
##

#######################
# The command interface!
#######################

library("optparse")
library(data.table)
library(lazyeval)
library(AnnotationDbi)

option_list = list(
  make_option(c("-p", "--platform"), type="character", default="", 
              help="Platform", metavar="character"),
  make_option(c("-i", "--inputFile"), type="character", default="", 
              help="inputFile", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

platform <- opt$platform
filePath <- opt$inputFile

# Load the platform
library(paste(platform, ".db", sep=""), character.only=TRUE)

# Read the data file
suppressWarnings(exprs <- fread(filePath, stringsAsFactors=FALSE, sep="\t", header=TRUE, autostart=10, data.table=FALSE, check.names=FALSE, fill=TRUE, na.strings="", showProgress=FALSE))
expr_ids <- exprs[colnames(exprs)[1]]

# Load these probes
x <- paste(platform, ".db", sep="")
probes <- AnnotationDbi::keys(get(x))

# Calculate the overlap
commonProbes <- intersect(unlist(expr_ids), probes)
percent <- ( length(commonProbes) / length(probes) ) * 100.0

# Send the result to stdout so parent process can pick it up
write(percent, stdout())
