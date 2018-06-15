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
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

platform <- opt$platform
filePath <- opt$inputFile
outFilePath <- opt$outputFile

# Get and use this DB
message("Loading...")
source("http://bioconductor.org/biocLite.R")

# library(Biobase)
# library("AnnotationDbi")

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

# find % overlap by column
overlap <- apply(index_data, 2, function(x) length(intersect(data$ID_REF, x)) / length(data$ID_REF))

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
colnames(data)[1] <- detected_id
#merged <- merge(data, index_data, by=detected_id, all.x=TRUE)

# Old version - R has 'choose the first' behavior
#out<-join(df1, df2, type="full")

# Take 2
# merged <- dplyr::inner_join(data, index_data, by=detected_id)
# message(head(merged))

converted_exprs <- index_data %>%  
  dplyr::select(c("ENSEMBL", detected_id)) %>%  # sub "ENSEMBL" for "GENEID"
  dplyr::inner_join(data, by = detected_id)

converted_exprs <- converted_exprs %>%
  dplyr::select(-rlang::UQ(rlang::sym(detected_id)))

message(head(converted_exprs))

# data <- merged[, c(3,2)]

# message("Replacing..")
# # Replace the probe IDs with Ensembles
# i <- 0
# c <- which( colnames(index_data) == detected_id )
# ec <- which( colnames(index_data) == "ENSEMBL" )
# for(m in data[, 1]){
# 	data[i, 1] = index_data[i, ec]
# 	i <- i + 1;
# }

# Remove all of the unmapped (NA) values (likely control probes?)
# data <- data[complete.cases(data), ]

# Data here can have duplicate rows that need to be squished together (via mean/max/etc),
# but this should be done at smash-time so that we can provide options on the squish method.

# Save to output file
write.table(converted_exprs, outFilePath, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
