# This script is based on the one written by J. Taroni:
# https://github.com/jaclyn-taroni/ref-txome/blob/master/4-athaliana_tximport.R
#
# Two required arguments:
# (1) --exp_dir    <experiement_dir> (directory where "quant.sf" files are found)
# (2) --gene2txmap <path_of_gene2txmap_file>
#
# Two Output files:
# (1) RDS file (saved as <experiment_dir>/txi_out.RDS)
# (2) TPM file (saved as <experiment_dir>/gene_lengthScaledTPM.tsv)

# Handle input arguments:
library("optparse")
option_list <- list(
  make_option("--exp_dir", type = "character"),
  make_option("--gene2txmap", type = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# directory of the experiment data
exp.directory <- opt$exp_dir

# Each sample's data is in a sub-dir of experiment directory.
samples <- list.dirs(exp.directory, full.names = FALSE, recursive = FALSE)

# Get a list of "quant.sf" files in experiment directory.
# The path of each "quant.sf" are in this format:
# <exp_dir>/<sample>/processed/quant.sf
sf.files <- list.files(exp.directory, pattern="^quant.sf$", full.names = TRUE, recursive = TRUE)

# Name each "sf.files" entry with the corresponding sample's name.
names(sf.files) <- samples

# Generate the gene to transcript mapping based on the input file
# specified by "--gene2txmap" option:
gene2tx.file <- opt$gene2txmap
gene2tx.df <- readr::read_tsv(gene2tx.file, col_names = c("gene_id", "tx_name", "DROP"))
gene2tx.df <- dplyr::select(gene2tx.df, -DROP)
# Reorder columns into: transcript, gene (which is required for tximport)
tx2gene <- gene2tx.df[, c("tx_name", "gene_id")]

# Run tximport and save RDS
txi.out <- tximport::tximport(files = sf.files,
                              type = "salmon",
                              txOut = TRUE,
                              tx2gene = tx2gene,
                              countsFromAbundance = "no")
rds_filename <- file.path(exp.directory, "txi_out.RDS")
saveRDS(txi.out, file = rds_filename)

# Run tximport::summarizeToGene() and save length-scaled
# TPM (transcripts per million) file.
txi.sum <- tximport::summarizeToGene(txi = txi.out,
                                     tx2gene = tx2gene,
                                     countsFromAbundance = "lengthScaledTPM")
# As data.frame, write to file
lstpm.df <- as.data.frame(txi.sum$counts)
lstpm.df <- cbind(rownames(lstpm.df), lstpm.df)
colnames(lstpm.df)[1] <- "Gene"
tpm_filename <- file.path(exp.directory, "gene_lengthScaledTPM.tsv")
readr::write_tsv(lstpm.df, path = tpm_filename)
