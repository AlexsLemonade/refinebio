# This script is based on the one written by J. Taroni:
# https://github.com/jaclyn-taroni/ref-txome/blob/master/4-athaliana_tximport.R
#
# Two required arguments:
# (1) --file_list    <path_of_file_list> (path containing list of quant.sf files)
# (2) --gene2txmap <path_of_gene2txmap_file>
# (3) --rds_file <path_of_rds_output>
# (3) --tpm_file <path_of_tmp_output>
#
# Two Output files:
# (1) RDS file (saved as value of --rds_file option)
# (2) Compressed TPM file (saved as --tpm_file option)

# Handle input arguments:
option_list <- list(
  optparse::make_option("--file_list", type = "character"),
  optparse::make_option("--gene2txmap", type = "character"),
  optparse::make_option("--rds_file", type = "character"),
  optparse::make_option("--tpm_file", type = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# text file that contains the full paths to quant.sf files
file_list <- opt$file_list
sf_files <- scan(file_list, character())

# We make the assumption that the quant paths are of a structure like:
# .../<SAMPLE_ACCESSION>/processed/quant.sf
samples <- lapply(sf_files, function(filename) {
    tokens <- unlist(strsplit(filename, "/")); tokens[length(tokens) - 1]
})
samples <- unlist(samples)

# Name each "sf_files" entry with the corresponding sample's name.
names(sf_files) <- samples

# Generate the gene to transcript mapping based on the input file
# specified by "--gene2txmap" option:
gene2tx_file <- opt$gene2txmap
gene2tx_df <- readr::read_tsv(gene2tx_file, col_names = c("gene_id", "tx_name"))
# Reorder columns into: transcript, gene (which is required for tximport)
tx2gene <- gene2tx_df[, c("tx_name", "gene_id")]

# Run tximport and save RDS
txi <- tximport::tximport(files = sf_files,
                          type = "salmon",
                          txOut = TRUE,
                          tx2gene = tx2gene,
                          countsFromAbundance = "no")
rds_filename <- opt$rds_file
saveRDS(txi, file = rds_filename)

# Run tximport::summarizeToGene() and save length-scaled
# TPM (transcripts per million) file.
txi_length_scaled <- tximport::tximport(files = sf_files,
                                        type = "salmon",
                                        tx2gene = tx2gene,
                                        countsFromAbundance = "lengthScaledTPM")

# As data.frame, write to file
lstpm_df <- as.data.frame(txi_length_scaled$counts)
lstpm_df <- tibble::rownames_to_column(lstpm_df, var = "Gene")
tpm_filename <- opt$tpm_file
readr::write_tsv(lstpm_df, path = tpm_filename)
