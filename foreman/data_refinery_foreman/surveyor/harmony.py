"""
Given a list of samples and their metadata, extract these common properties:

  `title`,
  `sex`,
  `age`,
  `specimen_part`,
  `genotype`,
  `disease`,
  `disease_stage`,
  `cell_line`,
  `treatment`,
  `race`,
  `subject`,
  `compound`,
  `time`

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
      '    ': 'donor B islets',
      'Technology Type': 'array assay',
      'Unit [time unit]': 'year',
      'sex': 'male'}

SRA Example:
    {"run_accession": "ERR1427135",
    "project_name": "Single-cell RNA sequencing of spleen-derived LCK cells from adult Zebrafish",
    "submission_accession": "ERA631093",
    "library_min_fragment_size": "",
    "bam_md5": "",
    "assembly_software": "",
    "library_prep_longitude": "",
    "library_selection": "cDNA",
    "pcr_isolation_protocol": "",
    "chip_protocol": "",
    "sequencing_primer_provider": "",
    "serotype": "",
    "environment_feature": "",
    "last_updated": "2018-11-16",
    "submitted_galaxy": "ftp.sra.ebi.ac.uk/vol1/run/ERR142/ERR1427135/LCK_4_90.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR142/ERR1427135/LCK_4_90.cram.crai",
    "extraction_protocol": "",
    "germline": "",
    "secondary_project": "",
    "culture_collection": "",
    "submission_tool": "",
    "sra_bytes": "",
    "read_strand": "",
    "rna_purity_280_ratio": "",
    "hi_c_protocol": "",
    "collected_by": "",
    "submitted_ftp": "ftp.sra.ebi.ac.uk/vol1/run/ERR142/ERR1427135/LCK_4_90.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR142/ERR1427135/LCK_4_90.cram.crai",
    "restriction_enzyme_target_sequence": "",
    "isolate": "",
    "fastq_bytes": "78258588;82464883",
    "instrument_platform": "ILLUMINA",
    "variety": "",
    "sequencing_date_format": "",
    "temperature": "",
    "sra_aspera": "",
    "ecotype": "",
    "submitted_aspera": "fasp.sra.ebi.ac.uk:/vol1/run/ERR142/ERR1427135/LCK_4_90.cram;fasp.sra.ebi.ac.uk:/vol1/run/ERR142/ERR1427135/LCK_4_90.cram.crai",
    "sampling_campaign": "",
    "bam_ftp": "",
    "tissue_lib": "",
    "environmental_sample": "",
    "control_experiment": "",
    "sex": "",
    "submitted_md5": "00c61590c4e662488b18a8397df876d6;ceb9b9368ba635b648d51bf15ad3626d",
    "checklist": "ERC000011",
    "fastq_galaxy": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/005/ERR1427135/ERR1427135_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/005/ERR1427135/ERR1427135_2.fastq.gz",
    "library_gen_protocol": "",
    "specimen_voucher": "",
    "library_prep_latitude": "",
    "submitted_bytes": "214700406;15844",
    "taxonomic_identity_marker": "",
    "run_date": "",
    "country": "",
    "ncbi_reporting_standard": "",
    "sample_description": "Protocols: The spleen from a heterozygote Tg(lck:EGFP) or wild-type fish was dissected and passed through a 40\u03bcm cell strainer using the plunger of a 1-mL syringe and cells were collected in cold 1xPBS/5% FBS. A non-transgenic line was used to set up the gating and exclude autofluorescent cells. Propidium iodide (PI) staining was used to exclude dead cells. Individual cells were sorted, using a Becton Dickinson Influx sorter with 488- and 561 nm lasers (Schulte et al., 2015) and collected in a single well of a 96 well plate containing 2.3 uL of 0.2 % Triton X-100 supplemented with 1 U/uL SUPERase In RNAse inhibitor (Ambion). The size, granularity and level of fluorescence for each cell were simultaneously recorded. The Smart-seq2 protocol (Picelli et al., 2014) was used to amplify the whole transcriptome and prepare libraries. Twenty-five cycles of PCR amplification were performed.",
    "experiment_title": "Illumina HiSeq 2500 paired end sequencing; Single-cell RNA sequencing of spleen-derived LCK cells from adult Zebrafish;Single-cell RNA sequencing of spleen-derived LCK cells from adult Zebrafish",
    "sra_galaxy": "",
    "sample_prep_interval": "",
    "fastq_md5": "d2aa4fe6ff43de680096df109a9e33f6;2664741d32d883f615205939ffc872d5",
    "sample_accession": "SAMEA4011859",
    "secondary_study_accession": "ERP015799",
    "experimental_protocol": "",
    "read_count": "4001164",
    "study_title": "Single-cell RNA sequencing of spleen-derived LCK cells from adult Zebrafish",
    "bio_material": "",
    "rna_prep_5_protocol": "",
    "host_body_site": "",
    "local_environmental_context": "",
    "assembly_quality": "",
    "collection_date_end": "",
    "sample_capture_status": "",
    "sample_title": "LCK_4#90",
    "host_genotype": "",
    "host_phenotype": "",
    "environmental_medium": "",
    "cultivar": "",
    "instrument_model": "Illumina HiSeq 2500",
    "faang_library_selection": "",
    "target_gene": "",
    "bam_bytes": "",
    "library_max_fragment_size": "",
    "experiment_target": "",
    "sequencing_date": "",
    "nominal_sdev": "20",
    "chip_ab_provider": "",
    "environment_material": "",
    "host_tax_id": "",
    "sample_material": "",
    "sample_storage_processing": "",
    "sra_md5": "",
    "cell_type": "",
    "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/005/ERR1427135/ERR1427135_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR142/005/ERR1427135/ERR1427135_2.fastq.gz",
    "sample_prep_interval_units": "",
    "broker_name": "ArrayExpress",
    "sub_strain": "",
    "base_count": "500145500",
    "library_strategy": "RNA-Seq",
    "restriction_site": "",
    "serovar": "",
    "investigation_type": "",
    "location": "",
    "library_source": "TRANSCRIPTOMIC",
    "sra_ftp": "",
    "library_layout": "PAIRED",
    "experimental_factor": "",
    "sequencing_primer_catalog": "",
    "environment_biome": "",
    "rna_purity_230_ratio": "",
    "dnase_protocol": "",
    "dev_stage": "",
    "library_prep_date_format": "",
    "bam_aspera": "",
    "binning_software": "",
    "rna_integrity_num": "",
    "library_prep_date": "",
    "location_start": "",
    "marine_region": "",
    "aligned": "",
    "file_location": "",
    "sample_collection": "",
    "chip_target": "",
    "nominal_length": "400",
    "broad_scale_environmental_context": "",
    "sequencing_location": "",
    "status": "public",
    "completeness_score": "",
    "lon": "",
    "fastq_aspera": "fasp.sra.ebi.ac.uk:/vol1/fastq/ERR142/005/ERR1427135/ERR1427135_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/ERR142/005/ERR1427135/ERR1427135_2.fastq.gz",
    "host_sex": "",
    "library_pcr_isolation_protocol": "",
    "sample_alias": "E-MTAB-4617:LCK_4#90",
    "mating_type": "",
    "collection_date_start": "",
    "sub_species": "",
    "contamination_score": "",
    "run_alias": "E-MTAB-4617:LCK_4#90",
    "restriction_enzyme": "",
    "depth": "",
    "submitted_read_type": "",
    "library_construction_protocol": "The spleen from a heterozygote Tg(lck:EGFP) or wild-type fish was dissected and passed through a 40\u03bcm cell strainer using the plunger of a 1-mL syringe and cells were collected in cold 1xPBS/5% FBS. A non-transgenic line was used to set up the gating and exclude autofluorescent cells. Propidium iodide (PI) staining was used to exclude dead cells. Individual cells were sorted, using a Becton Dickinson Influx sorter with 488- and 561 nm lasers (Schulte et al., 2015) and collected in a single well of a 96 well plate containing 2.3 uL of 0.2 % Triton X-100 supplemented with 1 U/uL SUPERase In RNAse inhibitor (Ambion). The size, granularity and level of fluorescence for each cell were simultaneously recorded.  The Smart-seq2 protocol (Picelli et al., 2014) was used to amplify the whole transcriptome and prepare libraries. Twenty-five cycles of PCR amplification were performed.",
    "host_growth_conditions": "",
    "collection_date": "",
    "experiment_alias": "E-MTAB-4617:LCK_4#90",
    "host_gravidity": "",
    "center_name": "Ludwig Center for Cancer Research, University of Lausanne, Lausanne, Switzerland  Swiss Institute of Bioinformatics",
    "identified_by": "",
    "cell_line": "",
    "sampling_site": "",
    "host": "",
    "library_name": "LCK_4#90",
    "tag": "ena;datahub",
    "first_created": "2016-05-26",
    "lat": "",
    "strain": "",
    "experiment_accession": "ERX1497667",
    "scientific_name": "Danio rerio",
    "host_status": "",
    "tax_id": "7955",
    "study_accession": "PRJEB14175",
    "submitted_format": "CRAM;CRAI",
    "submitted_host_sex": "",
    "bisulfite_protocol": "",
    "altitude": "",
    "rt_prep_protocol": "",
    "host_scientific_name": "",
    "bam_galaxy": "",
    "accession": "ERR1427135",
    "secondary_sample_accession": "ERS1182969",
    "sample_storage": "",
    "cage_protocol": "",
    "sampling_platform": "",
    "taxonomic_classification": "",
    "location_end": "",
    "protocol_label": "",
    "elevation": "",
    "salinity": "",
    "sequencing_method": "",
    "sequencing_primer_lot": "",
    "first_public": "2016-12-01",
    "transposase_protocol": "",
    "study_alias": "E-MTAB-4617",
    "library_prep_location": "",
    "rna_prep_3_protocol": "",
    "ph": "",
    "sequencing_longitude": "",
    "tissue_type": "Spleen",
    "isolation_source": "",
    "name": "E-MTAB-4617:LCK_4#90",
    "webinSubmissionAccountId": "Webin-24",
    "taxId": 7955,
    "release": "2016-12-01T17:01:24Z",
    "update": "2021-08-21T15:38:37.980Z",
    "submitted": "2016-05-26T17:59:08Z",
    "characteristics_ENA_first_public_0_text": "2016-12-01",
    "characteristics_ENA_first_public_0_tag": "attribute",
    "characteristics_ENA_last_update_0_text": "2016-05-26",
    "characteristics_ENA_last_update_0_tag": "attribute",
    "characteristics_ENA-CHECKLIST_0_text": "ERC000011",
    "characteristics_ENA-CHECKLIST_0_tag": "attribute",
    "characteristics_External_Id_0_text": "SAMEA4011859",
    "characteristics_External_Id_0_tag": "Namespace:BioSample",
    "characteristics_INSDC_center_alias_0_text": "Ludwig Center for Cancer Research, University of Lausanne, Lausanne, Switzerland Swiss Institute of Bioinformatics",
    "characteristics_INSDC_center_name_0_text": "Ludwig Center for Cancer Research, University of Lausanne, Lausanne, Switzerland Swiss Institute of Bioinformatics",
    "characteristics_INSDC_first_public_0_text": "2016-12-01T17:01:24Z",
    "characteristics_INSDC_last_update_0_text": "2016-05-26T17:59:08Z",
    "characteristics_INSDC_status_0_text": "public",
    "characteristics_SRA_accession_0_text": "ERS1182969",
    "characteristics_Submitter_Id_0_text": "E-MTAB-4617:LCK_4#90",
    "characteristics_Submitter_Id_0_tag": "Namespace:Ludwig Center for Cancer Research, University of Lausanne, Lausanne, Switzerland  Swiss Institute of Bioinformatics",
    "characteristics_broker_name_0_text": "ArrayExpress",
    "characteristics_broker_name_0_ontologyTerms_0": "http://www.ebi.ac.uk/efo/EFO_0009737",
    "characteristics_common_name_0_text": "zebrafish",
    "characteristics_description_0_text": "Protocols: The spleen from a heterozygote Tg(lck:EGFP) or wild-type fish was dissected and passed through a 40\u03bcm cell strainer using the plunger of a 1-mL syringe and cells were collected in cold 1xPBS/5% FBS. A non-transgenic line was used to set up the gating and exclude autofluorescent cells. Propidium iodide (PI) staining was used to exclude dead cells. Individual cells were sorted, using a Becton Dickinson Influx sorter with 488- and 561 nm lasers (Schulte et al., 2015) and collected in a single well of a 96 well plate containing 2.3 uL of 0.2 % Triton X-100 supplemented with 1 U/uL SUPERase In RNAse inhibitor (Ambion). The size, granularity and level of fluorescence for each cell were simultaneously recorded. The Smart-seq2 protocol (Picelli et al., 2014) was used to amplify the whole transcriptome and prepare libraries. Twenty-five cycles of PCR amplification were performed.",
    "characteristics_fsc_0_text": "23555",
    "characteristics_fsc_0_tag": "attribute",
    "characteristics_gfp_0_text": "109",
    "characteristics_gfp_0_tag": "attribute",
    "characteristics_individual_0_text": "2",
    "characteristics_individual_0_tag": "attribute",
    "characteristics_organism_0_text": "Danio rerio",
    "characteristics_organism_0_ontologyTerms_0": "http://purl.obolibrary.org/obo/NCBITaxon_7955",
    "characteristics_organism_1_text": "Danio rerio",
    "characteristics_organism_1_ontologyTerms_0": "http://purl.obolibrary.org/obo/NCBITaxon_7955",
    "characteristics_organism_1_tag": "attribute",
    "characteristics_pi_0_text": "9",
    "characteristics_pi_0_tag": "attribute",
    "characteristics_plate_0_text": "4",
    "characteristics_plate_0_tag": "attribute",
    "characteristics_ssc_0_text": "10994",
    "characteristics_ssc_0_tag": "attribute",
    "characteristics_tissue_0_text": "Spleen",
    "characteristics_tissue_0_ontologyTerms_0": "http://purl.obolibrary.org/obo/UBERON_0002106",
    "characteristics_tissue_0_tag": "attribute",
    "characteristics_title_0_text": "LCK_4#90",
    "characteristics_well_0_text": "B12",
    "characteristics_well_0_tag": "attribute",
    "relationships_0_source": "SAMEG321352",
    "relationships_0_type": "has member",
    "relationships_0_target": "SAMEA4011859",
    "externalReferences_0_url": "https://www.ebi.ac.uk/ena/browser/view/SAMEA4011859",
    "submittedVia": "JSON_API",
    "create": "2018-12-07T06:18:17.244Z",
    "_links_self_href": "https://www.ebi.ac.uk/biosamples/samples/SAMEA4011859",
    "_links_curationDomain_href": "https://www.ebi.ac.uk/biosamples/samples/SAMEA4011859{?curationdomain}",
    "_links_curationDomain_templated": true,
    "_links_curationLinks_href": "https://www.ebi.ac.uk/biosamples/samples/SAMEA4011859/curationlinks",
    "_links_curationLink_href": "https://www.ebi.ac.uk/biosamples/samples/SAMEA4011859/curationlinks/{hash}",
    "_links_curationLink_templated": true,
    "_links_structuredData_href": "https://www.ebi.ac.uk/biosamples/structureddata/SAMEA4011859",
    "study_ena_accession": "PRJEB14175",
    "study_ena_keywords": "",
    "study_ena_description": "Transcriptome data from individual lck:GFP expressing cells isolated from adult Zebrafish spleen. LCK is a marker of lymphocytes and here we identified two major subpopulations corresponding to T-cells and NK-like and a minor one of myeloid-like cells. Single cell transcriptomes are matched with FACS index sorting data (GFP, forward and side light scatter and dead cell staining)",
    "study_ena_study_name": "Single-cell RNA sequencing of spleen-derived LCK cells from adult Zebrafish",
    "study_ena_project_name": "Single-cell RNA sequencing of spleen-derived LCK cells from adult Zebrafish",
    "study_ena_geo_accession": "",
    "study_ena_isolate": "",
    "study_ena_center_name": "Ludwig Center for Cancer Research, University of Lausanne, Lausanne, Switzerland Swiss Institute of Bioinformatics",
    "study_ena_secondary_study_accession": "ERP015799",
    "study_ena_tag": "ena",
    "study_ena_secondary_study_center_name": "",
    "study_ena_strain": "",
    "study_ena_study_title": "Single-cell RNA sequencing of spleen-derived LCK cells from adult Zebrafish",
    "study_ena_last_updated": "2016-05-26",
    "study_ena_first_public": "2017-02-03",
    "study_ena_broker_name": "ArrayExpress",
    "study_ena_study_description": "Transcriptome data from individual lck:GFP expressing cells isolated from adult Zebrafish spleen. LCK is a marker of lymphocytes and here we identified two major subpopulations corresponding to T-cells and NK-like and a minor one of myeloid-like cells. Single cell transcriptomes are matched with FACS index sorting data (GFP, forward and side light scatter and dead cell staining)",
    "study_ena_secondary_study_alias": "E-MTAB-4617",
    "study_ena_scientific_name": "",
    "study_ena_study_alias": "E-MTAB-4617",
    "study_ena_breed": "",
    "study_ena_tax_id": "",
    "study_ena_parent_study_accession": "",
    "study_ena_study_accession": "PRJEB14175",
    "study_ena_submission_tool": "",
    "study_ena_cultivar": "",
    "study_ena_tax_division": "",
    "study_ena_status": "public"}

GEO:
    ex:
        {'channel_count': ['1'],
         'characteristics_ch1': ['patient: P-39',
                                 'gender: female',
                                 'age: 65',
                                 'location: lower leg',
                                 'transplanted organ: kidney',
                                 'immunosuppressive drugs: azathioprine + prednison',
                                 'sample type: squamous cell carcinoma',
                                 'cell type: keratinocyte'],
         'contact_address': ['Einthovenweg 20'],
         'contact_city': ['Leiden'],
         'contact_country': ['Netherlands'],
         'contact_department': ['Toxicogenetics, S4-P'],
         'contact_email': ['h.vrieling@lumc.nl'],
         'contact_institute': ['Leiden University Medical Center'],
         'contact_laboratory': ['room T4-34'],
         'contact_name': ['Harry,,Vrieling'],
         'contact_zip/postal_code': ['2333 ZC'],
         'data_processing': ['Raw data was extracted from the BeadChip data files in '
                             'Illumina’s BeadStudio Version 3.2 software using the '
                             'gene expression module (v 3.2.7). Background subtracted '
                             'data was further analyzed in R-based Bioconductor '
                             'package, lumi (version 1.12.4). In lumi, the data was '
                             'transformed (variance-stabilizing transformation (VST)) '
                             'and normalized (robust spline normalization (RSN)), '
                             'resulting in log-transformed normalized data. The '
                             'R-package illuminaHumanv2.db (version 1.4.1) was used '
                             'for annotation. The data were purged of genes that did '
                             'not meet the detection limit (expression-detection '
                             'P-value >0.01) and/or were not annotated. The limma '
                             'R-package (version 3.2.3) was used to identify '
                             'differentially expressed genes (DEGs) between SCC, AK '
                             'and NS. Gene set enrichment analysis (GSEA) was '
                             'performed with the significantly DEGs from the limma '
                             'analysis using DAVID Bioinformatic Resources v6.7 '
                             '(http://david.abcc.ncifcrf.gov). GSEA on the entire data '
                             'set was performed using the parametric gene set '
                             'enrichment analysis (PGSEA) R-package (version '
                             '1.14.0).      To identify activation of transcription '
                             'factors in AKs and SCCs, the DEGs from the limma '
                             'analysis were investigated using the online analysis '
                             'tool oPOSSUM.',
                             'Matrix normalized matrix shows VST-transformed, '
                             'RSN-normalized data (used scripts from lumi package)',
                             'Matrix non-normalized: AVG_Signal: average signal for '
                             'the probe; BEAD_STDERR: standard error of the beads; '
                             'Avg_NBEADS: average number of beads for that probe; '
                             'Detection Pval: detection p-value. All extracted from '
                             'Beadstudio'],
         'data_row_count': ['48701'],
         'description': ['1881436235_A'],
         'extract_protocol_ch1': ['RNA was isolated from SCC and AK samples that '
                                  'contained at least 70% tumor cells, as determined '
                                  'by haematoxylin and eosin stained frozen sections. '
                                  'From the sample of unexposed NS the epidermis was '
                                  'removed for further processing by cryosectioning '
                                  'parallel to the outer surface of the skin biopsy. '
                                  'RNA was extracted from frozen material using the '
                                  'RNeasy Fibrous Tissue kit (Qiagen), which included '
                                  'proteinase K treatment (10 min at 55˚C) of the '
                                  'lysed sample in RLT-buffer and on-column DNase '
                                  'treatment. RNA was quantified using a Nanodrop '
                                  '(NanoDrop technologies) and evaluated for '
                                  'degradation with a RNA 6000 Nano Labchip on the '
                                  '2100 Bioanalyzer (Agilent Technologies)'],
         'geo_accession': ['GSM808778'],
         'hyb_protocol': ['The standard Illumina hybridization protocol was used. In '
                          'brief, the samples were hybridized to the arrays at 58ºC '
                          'overnight.'],
         'label_ch1': ['biotin'],
         'label_protocol_ch1': ['100 ng of total RNA was converted to cDNA and '
                                'subsequently labeled cRNA using the Ambion Illumina '
                                'TotalPrep RNA Amplification kit (Ambion) according to '
                                'manufacturer’s instructions'],
         'last_update_date': ['Feb 06 2013'],
         'molecule_ch1': ['total RNA'],
         'organism_ch1': ['Homo sapiens'],
         'platform_id': ['GPL6102'],
         'scan_protocol': ['The beadChips were scanned using the Illumina BeadArray '
                           'Reader, using the standard Illumina scanning protocol'],
         'series_id': ['GSE32628', 'GSE32969', 'GSE32979'],
         'source_name_ch1': ['cutaneous squamous cell carcinoma'],
         'status': ['Public on Feb 06 2013'],
         'submission_date': ['Oct 05 2011'],
         'supplementary_file': ['NONE'],
         'taxid_ch1': ['9606'],
         'title': ['SCC_P-39'],
         'type': ['RNA']}
"""

import csv
import re
from io import StringIO
from typing import Dict, List

from data_refinery_common.logging import get_and_configure_logger
from data_refinery_foreman.surveyor.utils import requests_retry_session

logger = get_and_configure_logger(__name__)

# re.Pattern / re._pattern_type were not defined in 3.8?
RE_TYPE = re.compile("").__class__


TITLE_FIELDS = [
    "sample name",
    "sample title",
    "title",
    "subject number",
    "extract name",
    "labeled extract name",
]


def extract_title(sample: Dict, priority_field: str = None) -> str:
    """Given a flat sample dictionary, find the title"""
    if priority_field:
        title_fields = [priority_field] + [tf for tf in TITLE_FIELDS if tf != priority_field]
    else:
        title_fields = TITLE_FIELDS

    # Specifically look up for imported, non-SDRF AE samples
    for comment in sample.get("source_comment", []):
        if "title" in comment.get("name", ""):
            return comment["value"]

    expanded_title_fields = create_variants(title_fields)

    for title_field in expanded_title_fields:
        if title_field in sample:
            return sample[title_field]

    # If we can't even find a unique title for this sample
    # something has gone horribly wrong.
    return None


def create_variants(fields_list: List):
    """Given a list of strings, create variations likely to give metadata hits.

    Ex, given 'cell line', add the ability to hit on 'characteristic [cell_line]' as well.
    """
    variants = set()

    for field in fields_list:
        space_variants = [field]
        # Variate forms of multi-word strings
        if " " in field:
            space_variants.append(field.replace(" ", "_"))
            space_variants.append(field.replace(" ", "-"))
            space_variants.append(field.replace(" ", ""))

        for space_variant in space_variants:
            variants.add(space_variant)
            variants.add("characteristic [" + space_variant + "]")
            variants.add("characteristic[" + space_variant + "]")
            variants.add("characteristics [" + space_variant + "]")
            variants.add("characteristics[" + space_variant + "]")
            variants.add("comment [" + space_variant + "]")
            variants.add("comment[" + space_variant + "]")
            variants.add("comments [" + space_variant + "]")
            variants.add("comments[" + space_variant + "]")
            variants.add("factorvalue[" + space_variant + "]")
            variants.add("factor value[" + space_variant + "]")
            variants.add("factorvalue [" + space_variant + "]")
            variants.add("factor value [" + space_variant + "]")
            variants.add("sample_" + space_variant)
            variants.add("sample_host" + space_variant)
            variants.add("sample_sample_" + space_variant)  # Yes, seriously.

            # handle SRA flattened json response
            variants.add(re.compile("characteristics_" + space_variant + "(?:_\d+_text)?$"))

    return list(variants)


def parse_sdrf(sdrf_url: str) -> List:
    """Given a URL to an SDRF file, download parses it into JSON."""

    try:
        sdrf_response = requests_retry_session().get(sdrf_url, timeout=60)
    except Exception:
        logger.exception("Unable to fetch URL: " + sdrf_url)
        return []

    if sdrf_response.status_code != 200:
        logger.error("Unable to fetch URL: " + sdrf_url, response_code=sdrf_response.status_code)
        return []

    sdrf_text = sdrf_response.text

    samples = []

    reader = csv.reader(StringIO(sdrf_text), delimiter="\t")
    for offset, line in enumerate(reader):

        # Get the keys
        if offset == 0:
            keys = line
            continue

        sample_values = line

        # Skip malformed lines
        if len(sample_values) != len(keys):
            continue

        sample = {}
        for col, value in enumerate(sample_values):
            key = keys[col]
            sample[key.strip().lower()] = value
        samples.append(sample)

    return samples


def preprocess_geo_sample(sample) -> List:
    """
    Prepares items from GEO for harmonization
    """
    new_sample = {}
    for key, value in sample.metadata.items():

        if key == "characteristics_ch1":
            for pair in value:

                # This will almost always happen, except if we get
                # a malformed response from the server.
                if ":" in pair:
                    split = pair.split(":", 1)
                    new_sample[split[0].strip().lower()] = split[1].strip()
            continue

        # Probably won't be a list with length greater than one,
        # but maybe?
        new_sample[key.strip().lower()] = " ".join(value)

    return new_sample


class Harmonizer:
    def __init__(self):
        sex_fields = [
            "sex",
            "gender",
            "subject gender",
            "subject sex",
            # This looks reduntant, but there are some samples which use
            # Characteristic[Characteristic[sex]]
            "characteristic [sex]",
            "characteristics [sex]",
        ]
        self.sex_fields = create_variants(sex_fields)

        age_fields = [
            "age",
            "patient age",
            "age of patient",
            "age (years)",
            "age (yrs)",
            "age (months)",
            "age (days)",
            "age (hours)",
            "age at diagnosis",
            "age at diagnosis years",
            "age at diagnosis months",
            "age at diagnosis days",
            "age at diagnosis hours",
            "characteristic [age]",
            "characteristics [age]",
        ]
        self.age_fields = create_variants(age_fields)

        specimen_part_fields = [
            # AE
            "organism part",
            "cell type",
            "tissue",
            "tissue type",
            "tissue source",
            "tissue origin",
            "source tissue",
            "tissue subtype",
            "tissue/cell type",
            "tissue region",
            "tissue compartment",
            "tissues",
            "tissue of origin",
            "tissue-type",
            "tissue harvested",
            "cell/tissue type",
            "tissue subregion",
            "organ",
            "characteristic [organism part]",
            "characteristics [organism part]",
            # GEO
            "isolation source",
            "tissue sampled",
            "cell description",
        ]
        self.specimen_part_fields = create_variants(specimen_part_fields)

        genotype_fields = [
            "strain/background",
            "strain",
            "strain or line",
            "background strain",
            "genotype",
            "genetic background",
            "genetic information",
            "genotype/variation",
            "ecotype",
            "cultivar",
            "strain/genotype",
        ]
        self.genotype_fields = create_variants(genotype_fields)

        disease_fields = [
            "disease",
            "disease state",
            "disease status",
            "diagnosis",
            "infection with",
            "sample type",
        ]
        self.disease_fields = create_variants(disease_fields)

        disease_stage_fields = [
            "disease state",
            "disease staging",
            "disease stage",
            "grade",
            "tumor grade",
            "who grade",
            "histological grade",
            "tumor grading",
            "disease outcome",
            "subject status",
        ]
        self.disease_stage_fields = create_variants(disease_stage_fields)

        cell_line_fields = [
            "cell line",
        ]

        self.cell_line_fields = create_variants(cell_line_fields)

        treatment_fields = [
            "treatment",
            "treatment group",
            "treatment protocol",
            "drug treatment",
            "clinical treatment",
        ]
        self.treatment_fields = create_variants(treatment_fields)

        race_fields = [
            "race",
            "ethnicity",
            "race/ethnicity",
        ]
        self.race_fields = create_variants(race_fields)

        subject_fields = [
            # AE
            "subject",
            "subject id",
            "subject/sample source id",
            "subject identifier",
            "human subject anonymized id",
            "individual",
            "individual identifier",
            "individual id",
            "patient",
            "patient id",
            "patient identifier",
            "patient number",
            "patient no",
            "donor id",
            "donor",
            # SRA
            "source name",
        ]
        self.subject_fields = create_variants(subject_fields)

        developmental_stage_fields = [
            "developmental stage",
            "development stage",
            "development stages",
            # SRA.
            "dev stage",
        ]
        self.developmental_stage_fields = create_variants(developmental_stage_fields)

        compound_fields = [
            "compound",
            "compound1",
            "compound2",
            "compound name",
            "drug",
            "drugs",
            "immunosuppressive drugs",
        ]
        self.compound_fields = create_variants(compound_fields)

        time_fields = [
            "time",
            "initial time point",
            "start time",
            "stop time",
            "time point",
            "sampling time point",
            "sampling time",
            "time post infection",
        ]
        self.time_fields = create_variants(time_fields)

    def harmonize_value(self, field_name: str, value):
        if field_name == "age":
            try:
                return float(value)
            except ValueError:
                try:
                    return float(value.split(" ")[0])
                except ValueError:
                    # This is probably something weird, like a '.'
                    return
        elif field_name == "sex":
            if value.lower() in ["f", "female", "woman"]:
                return "female"
            elif value.lower() in ["m", "male", "man"]:
                return "male"
            else:
                return value.lower()
        else:
            return value.lower().strip()

    def harmonize_field(
        self, sample_metadata: Dict, harmonized_sample: Dict, variants_attribute: str
    ):
        sample_attribute = variants_attribute.split("_fields")[0]
        harmonized_values = set()
        key_variants = getattr(self, variants_attribute)
        key_regexes = [r for r in key_variants if isinstance(r, RE_TYPE)]
        key_strings = {s for s in key_variants if isinstance(s, str)}

        for key, value in sample_metadata.items():
            lower_key = key.lower().strip()
            if lower_key in key_strings:
                harmonized_value = self.harmonize_value(sample_attribute, value)
                if harmonized_value:
                    harmonized_values.add(harmonized_value)

            for key_regex in key_regexes:
                if key_regex.fullmatch(lower_key):
                    harmonized_value = self.harmonize_value(sample_attribute, value)
                    if harmonized_value:
                        harmonized_values.add(harmonized_value)

        # Age can only be one value in the database.
        if sample_attribute == "age":
            if len(harmonized_values) > 0:
                harmonized_sample[sample_attribute] = harmonized_values.pop()
        else:
            harmonized_sample[sample_attribute] = ";".join(
                sorted([str(v) for v in harmonized_values])
            )

    def harmonize_sample(self, sample_metadata: Dict, title_field: str = None) -> Dict:
        variants_attributes = [
            "sex_fields",
            "age_fields",
            "specimen_part_fields",
            "genotype_fields",
            "disease_fields",
            "disease_stage_fields",
            "cell_line_fields",
            "treatment_fields",
            "race_fields",
            "subject_fields",
            "developmental_stage_fields",
            "compound_fields",
            "time_fields",
        ]

        harmonized_sample = {}
        harmonized_sample["title"] = extract_title(sample_metadata, title_field)
        for variants_attribute in variants_attributes:
            self.harmonize_field(sample_metadata, harmonized_sample, variants_attribute)

        return harmonized_sample


def determine_title_field(a_samples: List[Dict], b_samples: List[Dict]) -> str:
    """Determines which field should be used for the title of the sample.

    Sometimes there is more metadata than actual samples, so we just
    take the field with the largest number of matching values.
    """
    max_title_field = ""
    max_matching_titles = 0

    # Reverse the order so if there's a tie for number of matches
    # we'll prioritize the fields at the front of the list.
    for title_field in reversed(TITLE_FIELDS):
        a_titles = {extract_title(sample, title_field) for sample in a_samples}
        b_titles = {extract_title(sample, title_field) for sample in b_samples}

        if a_titles == b_titles:
            max_title_field = title_field
            max_matching_titles = len(a_titles)

    if max_matching_titles > 0:
        return max_title_field


def harmonize_all_samples(sample_metadata: List[Dict], title_field: str = None) -> Dict:
    """Returns a mapping of sample title to harmonized sample metadata.

    See docstring at top of file for further clarfication of what "harmonized" means."""
    harmonizer = Harmonizer()

    harmonized_samples = {}
    for sample in sample_metadata:
        harmonized_sample = harmonizer.harmonize_sample(sample, title_field)
        title = harmonized_sample["title"]
        harmonized_samples[title] = harmonized_sample

    return harmonized_samples
