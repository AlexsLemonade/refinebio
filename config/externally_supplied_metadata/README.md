## Externally supplied metadata
#### E Flynn
#### Last updated 3/25/2020

All code to produce this update is included in the [sl_label repository](https://github.com/erflynn/sl_label repository)

This update includes imputed sex labels for microarray data (mouse, rat, and human).

The majority of gene expression data is missing metadata sex labels (see Table 1). This lack of labels prevents us from examining the breakdown by sex of many studies. We used the expression of X and Y chromosome genes and metadata sex labels to train a logistic regression model (with elastic net penalty) to predict sample sex. Across all three organisms, our models achieve approximately 95% accuracy in a randomly selected held-out test set as compared to the metadata labels. Additionally, we assessed the accuracy of our model, on various subsets of the data; comparing to all metadata sex labels (agreement 93.5-94.8%), a random sample of single sex studies (agreement 92.6-96.5%), and, in human, manually annotated sex labels from a previous analysis (94.2%) (see Table 2).  

| organism | Samples (n) | Samples missing metadata annotation | Studies (n) | Studies missing metadata annotation |
| ----- | ---- | ---- | ---- | ---- |
| human | 430119 | 74.90% | 14987 | 87.30% |
| mouse | 228707 | 74.00% | 12995 | 80.90% | 
| rat | 31361 | 55.50% | 1295 | 65.30% |
Table 1. Metadata missingness for sex labels.


| dataset | Human | Mouse | Rat |
| ----- | ---- | ---- | ---- | 
| Training (n=1400) | 95.60% | 96.10% | 98.60% |
| Testing (n=600) | 95.20% | 95.70% | 95.50% |
| Metadata | 93.5% (107748) | 93.5% (58473) | 94.8% (13995)|
| High-confidence | 97.3% (7301) | 95.5% (3968) | n/a |
| Single sex - f | 96.5% (12919) | 93.4% (13243) | 96.3% (2240) |
| Single sex - m | 92.6% (8128) | 93.6% (30225) | 95.8% (12689) |
| Manual annotations | 94.2% (8289) | n/a | n/a |

Table 2. Concordance of sex labels. Numbers in parentheses indicate the total number of samples, percentages the number of samples that agree divided by the total number of samples. High confidence labels have matching metadata and clustering based expression labels.


The cleaned metadata sex labels are also included in the `cleaned_metadata/` directory for microarray and RNA-seq. This process mapped all harmonized sex labels to "male", "female", "mixed", or "unknown". Code for this is included in the [`erflynn/sl_label`](https://github.com/erflynn/sl_label) repository under [`code/01_metadata/`](https://github.com/erflynn/sl_label/tree/bbc7f060a84598f48482f590c09b1df723d4d366/code/01_metadata). 
