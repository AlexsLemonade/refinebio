# refine.bio RNA-Seq Sample Compendium

This is a refine.bio RNA-seq sample compendium.
It represents the Salmon output for the collection of RNA-seq samples from an organism that have been processed with refine.bio.
Each individual sample has its own `quant.sf` file; the samples have not been aggregated and normalized.
RNA-seq sample compendia are designed to allow users that are comfortable handling these files to generate output that is most useful for their downstream applications.

Please see the [Salmon documentation on the `quant.sf` output format](https://salmon.readthedocs.io/en/latest/file_formats.html#quantification-file) for more information.

## Contents

This download includes individual `quant.sf` files that that are organized into to folders that represent experiments.
The individual files represent individual runs (e.g., SRR, ERR, DRR) and the folders represent projects (e.g., SRP, ERP, DRP) from Sequence Read Archive.
We also included any metadata available from refine.bio, which is limited for RNA-seq data at this time.

* The `aggregated_metadata.json` file contains experiment (project) metadata.

* Each folder corresponds to an experiment and contains a metadata TSV file and `quant.sf` files (named following the convention `<sample-accession-id>_quant.sf`) for each sample that refine.bio was able to process with Salmon.

* Sample metadata (e.g. disease vs. control labels) are contained in the TSV file with `metadata` in the filename as well as any JSON files.
We apply light harmonization to some sample metadata fields, which are denoted by `refinebio_` (`refinebio_annotations` is an exception).
The contents of a sample's `refinebio_annotations` field in JSON files include the submitter-supplied sample metadata.

* Experiment metadata (e.g., experiment title and description) are contained in JSON files with `metadata` in the filename.

_Note that we use the terms "sample" and "experiment" to be consistent with the rest of refine.bio, but files will use run identifiers and project identifiers, respectively._

Please see [our documentation](http://docs.refine.bio/en/latest/main_text.html#rna-seq-sample-compendia) for more details.

## Contact

If you identify issues with this download, please [file an issue on GitHub](https://github.com/AlexsLemonade/refinebio/issues).
If you would prefer to report issues via e-mail, you can also email [requests@ccdatalab.org](mailto:requests@ccdatalab.org).

## Citing refine.bio

Please use the following:

Casey S. Greene, Dongbo Hu, Richard W. W. Jones, Stephanie Liu, David S. Mejia, Rob Patro, Stephen R. Piccolo, Ariel Rodriguez Romero, Hirak Sarkar, Candace L. Savonen, Jaclyn N. Taroni, William E. Vauclain, Deepashree Venkatesh Prasad, Kurt G. Wheeler. **refine.bio: a resource of uniformly processed publicly available gene expression datasets.** URL: https://www.refine.bio

_Note that the contributor list is in alphabetical order as we prepare a manuscript for submission._
