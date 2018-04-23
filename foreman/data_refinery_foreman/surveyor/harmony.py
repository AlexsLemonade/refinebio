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
    
    Array Express Example:
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

    SRA Example:
        {'alias': 'GSM2997959_r1',
         'broker_name': 'GEO',
         'center_name': 'GEO',
         'center_project_name': 'GSE99065',
         'ena-base-count': '279111754922',
         'ena-spot-count': '1152605026',
         'experiment_accession': 'SRX3691797',
         'experiment_design_description': None,
         'experiment_title': 'NextSeq 500 paired end sequencing; GSM2997959: INOF_FRT; '
                             'Homo sapiens; RNA-Seq',
         'lab_name': '',
         'library_construction_protocol': 'cDNA library produced with TGIRT RNA was '
                                          'isolated and DNase treated using RNEasy '
                                          'mini kit (Qiagen 74106) according to the '
                                          'manufacturer protocol. 5ug DNA-free total '
                                          'RNA was then ribodepleted using Ribo-zero '
                                          'Gold (Illumina RZG1224 ) according to the '
                                          'manufacturer protocol and purified using a '
                                          'modified ZYMO RNA Clean and Concentrator '
                                          '(R1016) protocol where 8 volumes EtOH '
                                          'instead of 4. rRNA depleted RNA was '
                                          'fragmented with NEBNext Magnesium RNA '
                                          'Fragmentation Module (E6150) followed by '
                                          'dephosphorylation using T4PNK (mandel  ) '
                                          'and purified by same modified ZYMO '
                                          'protocol. cDNAs were synthesized via TGIRT '
                                          'template-switching with 1µM TGIRT-III '
                                          'reverse transcriptase (Ingex, LLC) for 15 '
                                          'min at 60o C, during which a DNA '
                                          'oligonucleotide containing the complement '
                                          'of an Illumina Read 2 sequencing '
                                          'primer-binding site becomes seamlessly '
                                          "linked to the 5' cDNA end. After reaction "
                                          'cleanup (Qiagen MinElute Reaction cleanup '
                                          "28206), a 5' adenylated DNA oligonucleotide "
                                          'containing the complement of an Illumina '
                                          'Read 1 sequencing primer-binding site is '
                                          "then ligated to the 3' cDNA end with "
                                          "Thermostable 5' AppDNA / RNA Ligase (New "
                                          'England Biolabs M0319). Properly ligated '
                                          'cDNAs were amplified by PCR (12 cycles) to '
                                          'synthesize the second strand and add '
                                          'Illumina flowcell capture and index '
                                          'sequences. Library was size-selected with '
                                          'Ampure XP beads (Beckman-Coulter) and '
                                          'quantified with Qubit and evaluated on an '
                                          'Agilent 2100 Bioanalyzer.',
         'library_layout': 'PAIRED',
         'library_selection': 'cDNA',
         'library_source': 'TRANSCRIPTOMIC',
         'library_strategy': 'RNA-Seq',
         'organism_id': '9606',
         'organism_name': 'HOMO SAPIENS',
         'platform_instrument_model': 'NextSeq500',
         'run_accession': 'SRR6718414',
         'run_ena_base_count': '1773379870',
         'run_ena_first_public': '2018-02-17',
         'run_ena_last_update': '2018-02-17',
         'run_ena_spot_count': '15046271',
         'sample_accession': 'SRS2951393',
         'sample_cell_type': 'Immortalized normal ovarian fibroblast',
         'sample_ena_base_count': '1773379870',
         'sample_ena_first_public': '2018-02-14',
         'sample_ena_last_update': '2018-02-14',
         'sample_ena_spot_count': '15046271',
         'sample_source_name': 'INOF cell line',
         'sample_title': 'INOF_FRT',
         'sample_treatment': 'none',
         'study_abstract': 'The ability to compare the abundance of one RNA molecule '
                           'to another is a crucial step for understanding how gene '
                           'expression is modulated to shape the transcriptome '
                           'landscape. However, little information is available about '
                           'the relative expression of the different classes of coding '
                           'and non-coding RNA or even between RNA of the same class. '
                           'In this study, we present a complete portrait of the human '
                           'transcriptome that depicts the relationship of all classes '
                           'of non-ribosomal RNA longer than sixty nucleotides. The '
                           'results show that the most abundant RNA in the human '
                           'rRNA-depleted transcriptome is tRNA followed by '
                           'spliceosomal RNA. Surprisingly, the signal recognition '
                           'particle RNA 7SL by itself occupied 8% of the ribodepleted '
                           'transcriptome producing a similar number of transcripts as '
                           'that produced by all snoRNA genes combined. In general, '
                           'the most abundant RNA are non-coding but many more protein '
                           'coding than non-coding genes produce more than 1 '
                           'transcript per million. Examination of gene functions '
                           'suggests that RNA abundance reflects both gene and cell '
                           'function. Together, the data indicate that the human '
                           'transcriptome is shaped by a small number of highly '
                           'expressed non-coding genes and a large number of '
                           'moderately expressed protein coding genes that reflect '
                           'cellular phenotypes. Overall design: RNA was isolated from '
                           'SKOV3ip1 and INOF human cell lines and selected with '
                           'different methods. The resulting libraries were '
                           'multiplexed and paired-end sequenced using Illumina HiSeq.',
         'study_accession': 'SRP107324',
         'study_ena_base_count': '279111754922',
         'study_ena_first_public': '2017-09-25',
         'study_ena_last_update': '2018-02-15',
         'study_ena_spot_count': '1152605026',
         'study_title': 'Simultaneous detection and relative quantification of coding '
                        'and non-coding RNA using a single sequencing reaction',
         'study_type': 'Transcriptome Analysis',
         'submission_accession': 'SRA562540',
         'submission_comment': 'submission brokered by GEO',
         'submission_title': 'Submitted by Gene Expression Omnibus on 25-SEP-2017'}

    """

    # Prepare the harmonized samples
    original_samples = metadata.copy()
    harmonized_samples = {}
    for sample in original_samples:
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip()
            if lower_key in age_fields:
                harmonized_samples[title]['age'] = int(value)
                break

    ##
    # Cell Parts!
    ##
    part_fields = [
                    # AE
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
                    'characteristics [organism part]',

                    # SRA
                    'sample_cell_type'
                ]
    part_fields = add_variants(part_fields)
    for sample in original_samples:
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
        for key, value in sample.copy().items():
            lower_key = key.lower().strip() 
            if lower_key in race_fields:
                harmonized_samples[title]['race'] = value

    ##
    # Subject
    ##
    subject_fields = [
                    # AE
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
                    'donor',

                    # SRA
                    'sample_source_name'
                ]
    subject_fields = add_variants(subject_fields)
    for sample in original_samples:
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
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
        title = sample['title']
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
