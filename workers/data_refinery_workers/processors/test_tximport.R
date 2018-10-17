option_list <- list(
  optparse::make_option("--txi_out", type = "character"),
  optparse::make_option("--gene2txmap", type = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

txi_out <- opt$txi_out
gene2txmap <- opt$gene2txmap

gene2tx <- readr::read_tsv(gene2txmap, col_names = c("gene_id", "tx_name"))
txi <- readRDS(txi_out)

# The rownames of the count matrix in txi will tell us about the whether or not
# the values are at the gene-level (what we want) or the transcript-level (what
# we do not want). In C. elegans, there's no overlap between the gene and tx
# identifiers.
gene_level <- any(gene2tx$gene_id %in% rownames(txi$counts))
transcript_level <- any(gene2tx$tx_name %in% rownames(txi$counts))

# if there are transcript ids in the rownames or no gene names in the rownames,
# throw an error
if (transcript_level | !(gene_level)) {
  stop("Count matrix is not at the gene level!")
}
