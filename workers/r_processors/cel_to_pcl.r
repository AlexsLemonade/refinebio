ProcessCelFiles <- function (process.dir, organism.code, output.file) {
  library("affy")
  library("affyio")
  # files now holds a list of all celfiles in the directory to be processed
  files <- list.celfiles(process.dir)


  # ptype is now a vector with the type of each array
  ptype <- sapply(files, function (f) read.celfile.header(paste(process.dir,
                                                                f,
                                                                sep="/"))[1])


  # ptype levels are the levels of that vector
  # (i.e. the types of arrays present) this is important because
  # only arrays of the same type can be processed together
  ptype.levels <- levels(as.factor(unlist(ptype)))


  # mapping holds the three major human platforms, the ReadAffy function
  # used below will read the dataset using the CDF specified in here.
  # This should correspond to the brainarray CDFs for Entrez.
  # I was only using three human platforms before
  # so they were hand defined but this loop should map them.
  mapping <- list()
  for (level in ptype.levels) {
    strip.version <- gsub('v2', '', level)
    strip.hyphen <- gsub('-', '', strip.version)
    stripped <- gsub('_', '', strip.hyphen)
    mapping[[level]] <- paste(stripped, organism.code, 'ENTREZG', sep='_')
  }


  for (level in ptype.levels) {
    # pfiles vector only contains files of this type
    pfiles <- paste(process.dir, subset(files, ptype==level), sep="/")
    # ReadAffy loads the array data using the custom CDF
    data <- ReadAffy(cdfname=mapping[[level]], filenames=pfiles)
    express <- rma(data)

    write.exprs(express, file=output.file)
  }
}
