processCelFiles <- function(processDir, organismCode, outputFile) {
    library("affy")
    library("affyio")
    ##files now holds a list of all celfiles in the directory to be processed
    files = list.celfiles(processDir)


    ##ptype is now a vector with the type of each array
    ptype = sapply(files, function(f) read.celfile.header(paste(processDir, f, sep="/"))[1])


    ##ptype levels are the levels of that vector (i.e. the types of arrays present)
    ##this is important because only arrays of the same type can be processed together
    ptype_levels = levels(as.factor(unlist(ptype)))


    ##mapping holds the three major human platforms, the ReadAffy function used below
    ##will read the dataset using the CDF specified in here.  This should correspond
    ##to the brainarray CDFs for Entrez.  I was only using three human platforms before
    ##so they were hand defined but this loop should map them.
    mapping <- list()
    for (level in ptype_levels) {
        stripversion = gsub('v2', '', level)
        striphyphen = gsub('-', '', stripversion)
        stripped = gsub('_', '', striphyphen)
        mapping[[level]] <- paste(stripped, organismCode, 'ENTREZG', sep='_')
    }


    for (level in ptype_levels) {
        ##pfiles vector only contains files of this type
        pfiles = paste(processDir, subset(files, ptype==level), sep="/")
        ##ReadAffy loads the array data using the custom CDF
        Data = ReadAffy(cdfname=mapping[[level]], filenames = pfiles)
        ##expresso processes the data
        ##express = expresso(Data, normalize.method="quantiles", bgcorrect.method="rma",pmcorrect.method="pmonly",summary.method="medianpolish")
        express = rma(Data)
        ##this line writes out the PCL file.  It will be named ENTREZG_$(level)_$(processDir).pcl
        ##write.exprs(express, file=paste(paste("ENTREZG", level, processDir, sep="_"), "pcl", sep="."))
        write.exprs(express, file=outputFile)
    }
}
