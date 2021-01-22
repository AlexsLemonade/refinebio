# refine.bio Normalized Compendium

This is a refine.bio normalized compendium comprised of all the samples from a species that we were able to process, aggregate, and normalize.
Normalized compendia provide a snapshot of the most complete collection of gene expression that refine.bio can produce for each supported organism.

You can read more about how we process refine.bio compendia in [our documentation](http://docs.refine.bio/en/latest/main_text.html#refine-bio-compendia).

## Contents

This download includes a gene expression matrix and experiment and sample metadata for all samples from a given organism that are fit for inclusion in the normalized compendium.

* The `aggregated_metadata.json` file contains experiment metadata and information about the transformation applied to the data.
Specifically, the `scale_by` field notes any row-wise transformation that was performed on the gene expression data. For normalized compendia, this value should always be `NONE`.

* The gene expression matrix is the tab-separated value (TSV) file that bears the species name.
For example, if you have downloaded the zebrafish normalized compendium, you would find the gene expression matrix in the file `DANIO_RERIO/DANIO_RERIO.tsv`.
Note that samples are _columns_ and rows are _genes_ or _features_.
This pattern is consistent with the input for many programs specifically designed for working with high-throughput gene expression data but may be transposed from what other machine learning libraries are expecting.

* Sample metadata (e.g. disease vs. control labels) are contained in the TSV file with `metadata` in the filename as well as any JSON files.
We apply light harmonization to some sample metadata fields, which are denoted by `refinebio_` (`refinebio_annotations` is an exception).
The contents of a sample's `refinebio_annotations` field include the submitter-supplied sample metadata.

* Experiment metadata (e.g., experiment title and description) are contained in JSON files with `metadata` in the filename.

Please see [our documentation](http://docs.refine.bio/en/latest/) for more details.

## Notes and observations

Combining all samples from a given species is a technical challenge, as it necessitates the integration of different microarray platforms and microarray data with RNA-seq data.
Although the normalization steps we perform eliminate some sources of technical bias, it is imperfect and an active area of development.
We strongly encourage you to consider using methods or models that can account for such biases and to explore and visualize the data with particular concern for samples' technology of origin (RNA-seq, microarray).

### Methods evaluation and exploratory data analysis

To identify appropriate methods for processing the initial releases of normalized compendia (described [here](http://docs.refine.bio/en/latest/main_text.html#species-compendia)), we performed a series of evaluations in a small zebrafish test compendium.
We've made these evaluations available and have documented our rationale on GitHub [here](https://github.com/AlexsLemonade/compendium-processing/tree/94089d2de170f0ca7b87e9e5c32239a8591faaa7/select_imputation_method).

We have also performed exploratory analyses in a larger zebrafish test compendium ([GitHub](https://github.com/AlexsLemonade/compendium-processing/tree/94089d2de170f0ca7b87e9e5c32239a8591faaa7/quality_check)).
We _briefly_ summarize our findings below, including links to relevant notebooks or plots:

* Genes that are longer tend to have higher values in RNA-seq data as compared to microarray data ([notebook](https://alexslemonade.github.io/compendium-processing/quality_check/07-technology_diff_exp.nb.html)).
* Unsurprisingly, shorter genes are less likely to be observed in RNA-seq data ([notebook](https://alexslemonade.github.io/compendium-processing/quality_check/06-lowly_expressed_genes.nb.html)).
* Genes that are often zero in RNA-seq data have lower average expression in microarray data ([notebook](https://alexslemonade.github.io/compendium-processing/quality_check/08-gene_lengths.nb.html)).
* We observe some differences in technology in the first two principle components, but there is also a group of RNA-seq samples that are different from all other samples (see below).
These are samples from the Wellcome Sanger Institute Zebrafish Mutation Project ([notebook](https://alexslemonade.github.io/compendium-processing/quality_check/11-rnaseq_bias.nb.html)).


![pca-test-compendium](https://raw.githubusercontent.com/AlexsLemonade/compendium-processing/6826cc448d8bd8605ba73d30e344e7d20438234c/quality_check/plots/larger_test_compendium_PCA.png)

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
