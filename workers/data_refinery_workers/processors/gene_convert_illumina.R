###
# gene_convert_illumina.R
# Alex's Lemonad Stand Foundation
# Childhood Cancer Data Lab
#
###

# import magrittr pipe
`%>%` <- dplyr::`%>%`

#######################
# The command interface!
#######################

library("optparse")
library(data.table)
library("dplyr")
library("rlang")
suppressPackageStartupMessages(library(AnnotationDbi))

option_list = list(
  make_option(c("-p", "--platform"), type="character", default="", 
              help="Platform", metavar="character"),
  make_option(c("-i", "--inputFile"), type="character", default="", 
              help="inputFile", metavar="character"),
  make_option(c("-o", "--outputFile"), type="character", default="", 
              help="outputFile", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

platform <- opt$platform
filePath <- opt$inputFile
outFilePath <- opt$outputFile

# Read the data file
message("Reading data file...")
suppressWarnings(exprs <- fread(filePath, 
					stringsAsFactors=FALSE, 
					sep="\t", header=TRUE, 
					autostart=10, 
					data.table=FALSE, 
					check.names=FALSE, 
					fill=TRUE, 
					na.strings="", 
					showProgress=FALSE)
				)

message("Here's the data..")
message(head(exprs))

# Build a PROBEID:ENSEMBL map
message("Here's the db..")
db_name <- paste(platform, ".db", sep="")
message(db_name)
library(db_name, character.only=TRUE)

# Replace the probe IDs with Ensembles
exprs[, 1] <- mapIds(get(db_name), keys=exprs[, 1], column="ENSEMBL", keytype="PROBEID")# , multiVals="first")
res <- exprs
message(head(res))

# Remove all of the unmapped (NA) values (likely control probes?)
res <- res[complete.cases(res), ]

# Data here can have duplicate rows that need to be squished together (via mean/max/etc),
# but this should be done at smash-time so that we can provide options on the squish method.

# Save to output file
write.table(res, outFilePath, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
