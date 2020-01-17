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

# need what the first column is called
exprs_id_name <- colnames(exprs)[1]
# don't replace the identifiers in exprs yet
mapped_list <- mapIds(get(db_name), keys=exprs[, 1], column="ENSEMBL", keytype="PROBEID", multiVals="list")
# get into data.frame form, should capture one to multiple mapping
mapped_df <- reshape2::melt(mapped_list)
# not 100% sure this bit is correct -- might be the other way around
colnames(mapped_df) <- c("ENSEMBL", "PROBEID")
annot_exprs_df <- mapped_df %>%
  dplyr::filter(!is.na(ENSEMBL)) %>%   # dropping unmapped -- mapIds might do this by default!
  dplyr::inner_join(y = exprs, by = c("PROBEID" = exprs_id_name)) %>%  # annotation heavy lifting
  dplyr::select(-PROBEID)  # drop the probe ids

# Remove named P Value columns:
if("Detection Pval" %in% colnames(annot_exprs_df))
{
	annot_exprs_df <- dplyr::select(annot_exprs_df, -c("Detection Pval"))
}
message(head(annot_exprs_df))

# Save to output file
write.table(annot_exprs_df, outFilePath, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
