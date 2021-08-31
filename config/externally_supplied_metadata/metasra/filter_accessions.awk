BEGIN {
    # The accessions file is tab-separated
    FS="\t"
    OFS="\t"

    print("Run Accession", "Sample Accession")
}

# Column #7 is "Type"
$7 == "RUN" {
    # Column #1 is the run's accession, and column #12 is the accession of the
    # associated sample
    print($1, $12)
}
