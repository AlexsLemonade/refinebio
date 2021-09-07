import csv
import json
import os
import sys
from typing import List

SRS_TO_SRR_MAP = None


def get_srr(srs_accession: str) -> List[str]:
    """MetaSRA uses the SRS accession for each sample, but we use the SRR on the backend"""

    global SRS_TO_SRR_MAP

    if SRS_TO_SRR_MAP is None:
        SRS_TO_SRR_MAP = {}
        with open(sys.argv[2], "r") as f:
            reader = csv.reader(f, delimiter="\t")
            # Skip the header row
            next(reader, None)
            for srr, srs in reader:
                try:
                    SRS_TO_SRR_MAP[srs].append(srr)
                except KeyError:
                    SRS_TO_SRR_MAP[srs] = [srr]

    out = SRS_TO_SRR_MAP.get(srs_accession, [])

    if len(out) == 0:
        print(f"No run accessions found for sample {srs_accession}")

    return out


def translate_attribute(metasra_attribute: dict) -> dict:
    """Translate a MetaSRA attribute that looks like this:
        {
            "property_id": "EFO:0000246",
            "unit_id": "missing",
            "value": 31.0
        }
    into our representation documented here:
    https://github.com/AlexsLemonade/refinebio/issues/2127#issuecomment-591651893
    """

    attribute = {"value": metasra_attribute["value"]}

    if metasra_attribute["unit_id"] != "missing":
        attribute["unit"] = metasra_attribute["unit_id"]

    return {metasra_attribute["property_id"]: attribute}


def translate_metasra_metadata(metasra_metadata: dict) -> (dict, dict):
    """Translate the MetaSRA json file into a dictionary of metadata and a
    dictionary of keywords that refine.bio can import
    """

    metadata = []
    keywords = {}

    counter = 0
    total = len(metasra_metadata)

    print(f"1/{total}")
    for srs_accession, value in metasra_metadata.items():
        counter += 1
        if counter % 5000 == 0:
            print(f"{counter}/{total}")

        accessions = get_srr(srs_accession)

        for accession in accessions:
            if len(value["real-value properties"]) != 0:
                metadata.append(
                    {
                        "sample_accession": accession,
                        "attributes": [
                            translate_attribute(a) for a in value["real-value properties"]
                        ],
                    }
                )

            if len(value["mapped ontology terms"]) != 0:
                keywords[accession] = value["mapped ontology terms"]

    print(f"{total}/{total}")
    return metadata, keywords


if __name__ == "__main__":
    with open(sys.argv[1], "r") as metasra:
        metasra_metadata = json.load(metasra)

    metadata, keywords = translate_metasra_metadata(metasra_metadata)

    with open("metasra_translated.json", "w+") as output:
        json.dump(metadata, output, indent=4)

    with open("metasra_keywords.json", "w+") as output:
        json.dump(keywords, output, indent=4)
