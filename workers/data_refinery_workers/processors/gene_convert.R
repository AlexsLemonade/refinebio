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
  make_option(c("-p", "--platform"), type="character", default="", 
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

library(Biobase)
library("AnnotationDbi")

# Read the data file
message("Reading data file...")
suppressWarnings(data <- fread(filePath, 
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
index_path = paste0('zcat /home/user/data_store/gene_indexes/', platform, '.tsv.gz')
suppressWarnings(index_data <- fread(index_path, 
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

# we'll call these the supported identifiers
supported_ids <- c("PROBEID", "ENTREZID", "SYMBOL", "ENSEMBL", "UNIGENE")

# message(index_data)
# res <- select(index_data, keys=keys(index_data), columns=supported_ids)

# find % overlap by column
overlap <- apply(index_data, 2, function(x) length(intersect(data$ID_REF, x)) / length(data$ID_REF))

# threshold here can be adjusted
# if (any(overlap > 0.9)) {
if (any(overlap > 0.6)) {
	# take column name with highest overlap value 
  	# overlap is a named numeric vector
  	detected_id <- names(which.max(overlap))
} else {
  	stop("Not enough overlapping ids detected!")
}

# Replace the probe IDs with Ensembles
i <- 0
c <- which( colnames(index_data) == detected_id )

for(m in data[, 1]){
	data[i, 1] = index_data[i, c]
	i <- i + 1;
}

# Remove all of the unmapped (NA) values (likely control probes?)
data <- data[complete.cases(data), ]

# Data here can have duplicate rows that need to be squished together (via mean/max/etc),
# but this should be done at smash-time so that we can provide options on the squish method.

# Save to output file
write.table(data, outFilePath, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
