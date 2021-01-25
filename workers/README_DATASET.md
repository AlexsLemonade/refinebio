# refine.bio Aggregated Dataset

This is a refine.bio dataset.

## Contents

This download includes gene expression matrices and experiment and sample metadata for the samples that you selected for download.

* The `aggregated_metadata.json` file contains information about the options you selected for download.
Specifically, the `aggregate_by` and `scale_by` fields note how the samples are grouped into gene expression matrices and how the gene expression data values were transformed, respectively.
The `quantile_normalized` fields notes whether or not quantile normalization was performed.
Currently, we only support skipping quantile normalization for RNA-seq experiments when aggregating by experiment on the web interface.

* Individual gene expression matrices and their corresponding sample metadata files are in their own directories.

* Gene expression matrices are the tab-separated value (TSV) files named by the experiment accession number (if aggregated by experiment) or species name (if aggregated by species).
Note that samples are _columns_ and rows are _genes_ or _features_.
This pattern is consistent with the input for many programs specifically designed for working with high-throughput gene expression data but may be transposed from what other machine learning libraries are expecting.

* Sample metadata (e.g. disease vs. control labels) are contained in TSV files with `metadata` in the filename as well as any JSON files.
We apply light harmonization to some sample metadata fields, which are denoted by `refinebio_` (`refinebio_annotations` is an exception).
The contents of a sample's `refinebio_annotations` field include the submitter-supplied sample metadata.

* Experiment metadata (e.g., experiment title and description) are contained in JSON files with `metadata` in the filename.

Please see [our documentation](https://refinebio-docs.readthedocs.io/) for more details.

## Usage

The gene expression matrix TSV and JSON files can be read in, manipulated, or parsed with standard functions or libraries in the language of your choice.
Below are some code snippets to help you import the data into R or Python and examine it.

### Reading TSV Files

Here's an example reading a gene expression TSV (`GSE11111.tsv`) into R as a data.frame with base R:

```
expression_df <- read.delim("GSE11111.tsv", header = TRUE,
							row.names = 1, stringsAsFactors = FALSE)
```

### Reading JSON Files

#### R

The `rjson` R package allows us to read a metadata JSON file (`aggregated_metadata.json`) into R as a list:

```
library(rjson)
metadata_list <- fromJSON(file = "aggregated_metadata.json")
```

#### Python

In Python, we can read in the metadata JSON like so:

```
import json
with open('aggregated_metadata.json', 'r') as jsonfile:
    data = json.load(jsonfile)
print(data)
```

For example R workflows, such as clustering of gene expression data, please see [our repository of example uses](https://github.com/AlexsLemonade/refinebio-examples).

## Contact

If you identify issues with this download, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues).
If you would prefer to report issues via e-mail, you can also email [requests@ccdatalab.org](mailto:requests@ccdatalab.org).

## Citing refine.bio

Please use the following:

Casey S. Greene, Dongbo Hu, Richard W. W. Jones, Stephanie Liu, David S. Mejia, Rob Patro, Stephen R. Piccolo, Ariel Rodriguez Romero, Hirak Sarkar, Candace L. Savonen, Jaclyn N. Taroni, William E. Vauclain, Deepashree Venkatesh Prasad, Kurt G. Wheeler. **refine.bio: a resource of uniformly processed publicly available gene expression datasets.** URL: https://www.refine.bio

_Note that the contributor list is in alphabetical order as we prepare a manuscript for submission._
