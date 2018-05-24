# This script is based on the one written by J. Taroni:
# https://github.com/jaclyn-taroni/ref-txome/blob/master/4-athaliana_tximport.R
#
# Two required arguments:
# (1) --exp_dir    <experiement_dir> (directory where "quant.sf" files are found)
# (2) --gene2txmap <path_of_gene2txmap_file>
#
# Two Output files:
# (1) RDS file (saved as <experiment_dir>/txi_out.RDS)
# (2) Compressed TPM file (saved as <experiment_dir>/gene_lengthScaledTPM.tsv.gz)

# Handle input arguments:
option_list <- list(
  optparse::make_option("--exp_dir", type = "character"),
  optparse::make_option("--gene2txmap", type = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# directory of the experiment data
exp_dir <- opt$exp_dir

# Each sample's data is in a sub-dir of experiment directory.
samples <- list.dirs(exp_dir, full.names = FALSE, recursive = FALSE)

# Get a list of "quant.sf" files in experiment directory.
# The paths of "quant.sf" files are in this format:
# <exp_dir>/<sample>/processed/quant.sf
sf_files <- list.files(exp_dir, pattern="^quant.sf$", full.names = TRUE, recursive = TRUE)

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
                          txOut = FALSE,
                          tx2gene = tx2gene,
                          countsFromAbundance = "no")
rds_filename <- file.path(exp_dir, "txi_out.RDS")
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
tpm_filename <- file.path(exp_dir, "gene_lengthScaledTPM.tsv")
readr::write_tsv(lstpm_df, path = tpm_filename)
# Compress this output file
R.utils::gzip(tpm_filename)
