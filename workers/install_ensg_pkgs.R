install.packages("xml2")
library("xml2")
ensg_url <- "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/22.0.0/ensg.asp"
html_content <- read_html(ensg_url)

# The second table in the html is the main table that lists all ensg packages.
data_table <- xml_find_all(html_content, ".//table")[[2]]

# Extract data rows from the second table:
data_rows <- xml_children(data_table)[3:xml_length(data_table)]

# These two global variables will be modified by save_chip_pkg():
chips <- c()
pkg_urls <- c()

# This function parses a data row in the table and saves chip names and
# URLs of "P" R source package to chips and pkg_urls respectively.
save_chip_pkg <- function(row) {
    # Column #3: "Chip"
    chip_name <- xml_text(xml_child(row, 3))
    # Convert chip_name to lower case, then normalize it by removing
    # "-", "_" and " " characters.
    chip_name <- gsub("[-_ ]", "", tolower(chip_name))

    # Column #16: "R Source Packages (C and P)".
    # We need the URL of "P", which is the second URL.
    source_pkgs <- xml_child(row, 16)
    probe_pkg_url <- xml_attr(xml_child(source_pkgs, 2), "href")

    if (nzchar(chip_name) && nzchar(probe_pkg_url)) {
        chips <<- c(chips, chip_name)
        pkg_urls <<- c(pkg_urls, probe_pkg_url)
    }
}

# Extract data of interest out of each row
lapply(data_rows, save_chip_pkg)

# Write chips and pkg_urls to a tab-delimited file
output_filename <- "/home/user/r_ensg_probe_pkgs.txt"
write.table(list(chips, pkg_urls), file=output_filename, quote=FALSE,
            row.names=FALSE, col.names=FALSE, sep="\t")

# Install these ensg packages
lapply(pkg_urls, devtools::install_url)
