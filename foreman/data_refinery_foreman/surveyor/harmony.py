from data_refinery_common.models import (
    Sample)
import requests

def add_variants(original_list):
    """ Adds variants to a list """
    copy = original_list.copy()
    for item in original_list:
        if ' ' in item:
            copy.append(item.replace(' ', '_'))
            copy.append(item.replace(' ', '-'))
            copy.append(item.replace(' ', ''))
        copy.append("characteristic [" + item + "]")
        copy.append("characteristic[" + item + "]")
        copy.append("characteristics [" + item + "]")
        copy.append("characteristics[" + item + "]")
        copy.append("comment [" + item + "]")
        copy.append("comment[" + item + "]")
        copy.append("comments [" + item + "]")
        copy.append("comments[" + item + "]")
        copy.append("factorvalue[" + item + "]")
        copy.append("factor value[" + item + "]")
        copy.append("factorvalue [" + item + "]")
        copy.append("factor value [" + item + "]")
    return copy

def harmonize(metadata):
    """ 

    Given some samples and metadata, harmonize into something universal.
    

 {'Array Data File': 'C30061.CEL',
  'Array Design REF': 'A-AFFY-1',
  'Assay Name': '1009003-C30061',
  'Characteristics[age]': '38',
  'Characteristics[developmental stage]': 'adult',
  'Characteristics[organism part]': 'islet',
  'Characteristics[organism]': 'Homo sapiens',
  'Characteristics[sex]': 'male',
  'Comment [ArrayExpress FTP file]': 'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.raw.1.zip',
  'Comment [Derived ArrayExpress FTP file]': 'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3050/E-MTAB-3050.processed.1.zip',
  'Derived Array Data File': 'C30061.txt',
  'Description': 'Islets from 38 years old male. Time from islet preparation '
                 'to culture initiation: 96 hours. Provided by the North West '
                 'Tissue Center Seattle.',
  'Extract Name': 'donor B differentiated cells RNA',
  'Factor Value[cell type]': 'differentiated',
  'Factor Value[individual]': 'B',
  'Factor Value[test result]': 'unsuccessful',
  'Image File': 'C30061.DAT',
  'Label': 'biotin',
  'Labeled Extract Name': 'donor B differentiated cells LEX',
  'Material Type': 'cell',
  'Protocol REF': 'P-MTAB-41862',
  'Source Name': 'donor B islets',
  'Technology Type': 'array assay',
  'Unit [time unit]': 'year',
  'sex': 'male'}

    """

    # Prepare the harmonized samples
    original_samples = metadata.copy()
    harmonized_samples = {}
    for sample in original_samples:
        title = sample['Assay Name']
        harmonized_samples[title] = {}

    ##
    # Sex!
    ##
    sex_fields = [    
                    'sex',
                    'gender',
                    'subject gender',
                    'subjext sex',
                    'characteristic [sex]',
                    'characteristics [sex]',
                   ]
    sex_fields = add_variants(sex_fields)

    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip()
            if lower_key in sex_fields:
                if "f" in value:
                    harmonized_samples[title]['sex'] = "female"
                    break
                else:
                    harmonized_samples[title]['sex'] = "male"
                    break

    ##
    # Age!
    ##
    age_fields = [    
                    'age',
                    'patient age',
                    'age of patient',
                    'age (years)',
                    'age at diagnosis',
                    'age at diagnosis years',
                    'characteristic [age]',
                    'characteristics [age]',
                ]
    age_fields = add_variants(age_fields)

    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip()
            if lower_key in age_fields:
                harmonized_samples[title]['age'] = int(value)
                break

    ##
    # Cell Parts!
    ##
    part_fields = [
                    'organism part', 
                    'cell type', 
                    'tissue', 
                    'tissue type', 
                    'tissue source', 
                    'tissue origin', 
                    'source tissue', 
                    'tissue subtype', 
                    'tissue/cell type', 
                    'tissue region', 
                    'tissue compartment', 
                    'tissues', 
                    'tissue of origin', 
                    'tissue-type', 
                    'tissue harvested', 
                    'cell/tissue type', 
                    'tissue subregion', 
                    'organ',
                    'characteristic [organism part]',
                    'characteristics [organism part]'
                ]
    part_fields = add_variants(part_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in part_fields:
                harmonized_samples[title]['part'] = value
                break

    ##
    # Genotype!
    ##
    genotype_fields = [
                    'strain/background', 
                    'strain', 
                    'strain or line', 
                    'background strain', 
                    'genotype', 
                    'genetic background', 
                    'genotype/variation', 
                    'ecotype', 
                    'cultivar', 
                    'strain/genotype'
                ]
    genotype_fields = add_variants(genotype_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in genotype_fields:
                harmonized_samples[title]['genotype'] = value

    ##
    # Disease!
    ##
    disease_fields = [
                    'disease', 
                    'disease state', 
                    'disease status', 
                    'diagnosis', 
                ]
    disease_fields = add_variants(disease_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in disease_fields:
                harmonized_samples[title]['disease'] = value

    disease_stage_fields = [
                    'disease staging', 
                    'disease stage', 
                    'grade', 
                    'tumor grade', 
                    'who grade', 
                    'histological grade', 
                    'tumor grading'
                ]
    disease_stage_fields = add_variants(disease_stage_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in disease_stage_fields:
                harmonized_samples[title]['disease_stage'] = value

    ##
    # Cell Line!
    ##
    cell_line_fields = [
                    'cell_line'
                ]
    cell_line_fields = add_variants(cell_line_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in cell_line_fields:
                harmonized_samples[title]['cell_line'] = value

    ##
    # Treatment!
    ##
    treatment_fields = [
                    'treatment', 
                    'treatment group', 
                    'treatment protocol', 
                    'drug treatment', 
                    'clinical treatment',
                ]
    treatment_fields = add_variants(treatment_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in treatment_fields:
                harmonized_samples[title]['treatment'] = value
    ##
    # Race!
    ##
    race_fields = [
                    'race', 
                    'ethnicity', 
                    'race/ethnicity'
                ]
    race_fields = add_variants(race_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in race_fields:
                harmonized_samples[title]['race'] = value

    ##
    # Subject
    ##
    subject_fields = [
                    'subject', 
                    'subject id', 
                    'subject/sample source id', 
                    'subject identifier', 
                    'human subject anonymized id', 
                    'individual', 
                    'individual identifier', 
                    'individual id', 
                    'patient', 
                    'patient id', 
                    'patient identifier', 
                    'patient number', 
                    'patient no', 
                    'donor id', 
                    'donor'
                ]
    subject_fields = add_variants(subject_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in subject_fields:
                harmonized_samples[title]['subject'] = value

    ##
    # Developement Stage!
    ##
    development_stage_fields = [
                    'developmental stage', 
                    'development stage', 
                    'development stages'
                ]
    development_stage_fields = add_variants(development_stage_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in development_stage_fields:
                harmonized_samples[title]['developmental_stage'] = value

    ##
    # Compound!
    ##
    compound_fields = [
                    'compound', 
                    'compound1', 
                    'compound2', 
                    'compound name', 
                    'drug', 
                    'drugs'
                ]
    compound_fields = add_variants(compound_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in compound_fields:
                harmonized_samples[title]['compound'] = value

    ##
    # Time!
    ##
    time_fields = [
                    'time', 
                    'initial time point', 
                    'start time', 
                    'stop time', 
                    'time point', 
                    'sampling time point', 
                    'sampling time'
                ]
    time_fields = add_variants(time_fields)
    for sample in original_samples:
        title = sample['Assay Name']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in time_fields:
                harmonized_samples[title]['time'] = value

    return harmonized_samples


def parse_sdrf(sdrf_url):
    """ Given a URL to an SDRF file, parsers it into JSON """

    try:
        sdrf_text = requests.get(sdrf_url, timeout=5).text
    except Exception:
        return []

    samples = []

    lines = sdrf_text.split('\n')
    for offset, line in enumerate(lines):

        # Get the keys
        if offset == 0:
            keys = line.split('\t')
            continue

        # Skip blank lines
        if line == "":
            continue

        sample_values = line.split('\t')
        if len(sample_values) != len(keys):
            import pdb
            pdb.set_trace()
            continue

        sample = {}
        for col, value in enumerate(sample_values):
            key = keys[col]
            sample[key] = value
        samples.append(sample)

    return samples
