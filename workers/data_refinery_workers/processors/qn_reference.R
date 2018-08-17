###
# qn_reference
# Alex's Lemonad Stand Foundation
# Childhood Cancer Data Lab
#
###

#######################
# The command interface!
#######################

library("optparse")
library(data.table)
library("preprocessCore")

option_list = list(
  make_option(c("-i", "--inputFile"), type="character", default="", 
              help="inputFile", metavar="character"),
  make_option(c("-o", "--outputFile"), type="character", default="", 
              help="outputFile", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

filePath <- opt$inputFile
outFilePath <- opt$outputFile

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

# Create the target result
data_matrix <- data.matrix(data, rownames.force = NA)
data_matrix<-data_matrix[,-1] # We don't want the index column

target <- normalize.quantiles.determine.target(data_matrix)

# Save to output file
message("Writing target!")
write.table(target, outFilePath, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
