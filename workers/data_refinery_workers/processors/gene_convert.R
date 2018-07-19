###
# gene_convert.R
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

# Read the data file
message("Reading master index...")
index_path = paste0('zcat /home/user/gene_indexes/', platform, '.tsv.gz')
suppressWarnings(index_exprs <- fread(index_path, 
					stringsAsFactors=FALSE, 
					sep="\t", 
					header=TRUE, 
					autostart=10, 
					data.table=FALSE, 
					check.names=FALSE, 
					fill=TRUE, 
					na.strings="", 
					showProgress=FALSE)
				)

# find % overlap by column
overlap <- apply(index_exprs, 2, function(x) length(intersect(exprs$ID_REF, x)) / length(exprs$ID_REF))

# threshold here can be adjusted
if (any(overlap >= 0.5)) {
	# take column name with highest overlap value 
  	# overlap is a named numeric vector
  	detected_id <- names(which.max(overlap))
} else {
	# Less than 50%? Get outta here!
  	stop("Not enough overlapping ids detected!")
}

message("Merging..")
converted_exprs <- index_exprs %>%
    dplyr::select(c("ENSEMBL", detected_id)) %>%
    dplyr::inner_join(exprs, by = detected_id) %>%
    dplyr::select(-rlang::UQ(rlang::sym(detected_id)))

# Data here can have duplicate rows that need to be squished together (via mean/max/etc),
# but this should be done at smash-time so that we can provide options on the squish method.

# Save to output file
write.table(converted_exprs, outFilePath, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
