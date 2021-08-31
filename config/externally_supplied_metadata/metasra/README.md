# Generating new metadata

To generate new metadata, it should be as simple as changing the URL in `translate.sh` to point to the new version of the MetaSRA JSON file and then re-running `translate.sh`. After this is done, upload `metasra_translated.json` and `metasra_keywords.json` to the `data-refinery-test-assets` S3 bucket so that we can use them to ingest new metadata.

# Importing metadata

Once you have uploaded the new metadata, you can import it with the following commands on the foreman instance:

First, you need to import the relevant ontologies if you haven't already:

```sh
for ontology in CL CVCL DOID EFO UBERON UO; do ./run_management_command.sh import_ontology --ontology "$ontology"; done
```

then you can import the MetaSRA metadata:

```sh
./run_management_command.sh import_external_sample_attributes --source-name "MetaSRA" --methods-url "https://pubmed.ncbi.nlm.nih.gov/28535296/" --file "s3://data-refinery-test-assets/metasra_translated.json"

./run_management_command.sh import_external_sample_keywords --source-name "MetaSRA" --methods-url "https://pubmed.ncbi.nlm.nih.gov/28535296/" --file "s3://data-refinery-test-assets/metasra_keywords.json"
```

Note that it is safe to run both of these commands again when there is new metadata to ingest, as long as you make sure that the source name and methods URL are identical.
