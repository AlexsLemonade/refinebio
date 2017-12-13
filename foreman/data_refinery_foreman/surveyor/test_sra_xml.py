EXPERIMENT_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="DRX001563&amp;display=xml">
<EXPERIMENT accession="DRX001563" center_name="RIKEN_CDB" alias="DRX001563" broker_name="DDBJ">
     <IDENTIFIERS>
          <PRIMARY_ID>DRX001563</PRIMARY_ID>
          <SUBMITTER_ID namespace="RIKEN_CDB">DRX001563</SUBMITTER_ID>
     </IDENTIFIERS>
     <TITLE>Illumina HiSeq 2000 sequencing; Exp_Gg_HH16_1_embryo_mRNAseq</TITLE>
     <STUDY_REF accession="DRP000595" refcenter="RIKEN_CDB" refname="DRP000595">
          <IDENTIFIERS>
               <PRIMARY_ID>DRP000595</PRIMARY_ID>
               <EXTERNAL_ID namespace="BioProject" label="BioProject ID">PRJDB90</EXTERNAL_ID>
               <SUBMITTER_ID namespace="RIKEN_CDB">DRP000595</SUBMITTER_ID>
          </IDENTIFIERS>
     </STUDY_REF>
     <DESIGN>
          <DESIGN_DESCRIPTION>Experiment for mRNAseq of chicken at stage HH16 (biological replicate 1)</DESIGN_DESCRIPTION>
          <SAMPLE_DESCRIPTOR accession="DRS001521" refcenter="RIKEN_CDB" refname="DRS001521">
               <IDENTIFIERS>
                    <PRIMARY_ID>DRS001521</PRIMARY_ID>
                    <EXTERNAL_ID namespace="BioSample" label="BioSample ID">SAMD00004104</EXTERNAL_ID>
                    <SUBMITTER_ID namespace="RIKEN_CDB">DRS001521</SUBMITTER_ID>
               </IDENTIFIERS>
          </SAMPLE_DESCRIPTOR>
          <LIBRARY_DESCRIPTOR>
               <LIBRARY_NAME>Gg_HH16_1_embryo_mRNAseq</LIBRARY_NAME>
               <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
               <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
               <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
               <LIBRARY_LAYOUT>
                    <SINGLE/>
               </LIBRARY_LAYOUT>
          </LIBRARY_DESCRIPTOR>
          <SPOT_DESCRIPTOR>
               <SPOT_DECODE_SPEC>
                    <SPOT_LENGTH>100</SPOT_LENGTH>
                    <READ_SPEC>
                         <READ_INDEX>0</READ_INDEX>
                         <READ_CLASS>Application Read</READ_CLASS>
                         <READ_TYPE>Forward</READ_TYPE>
                         <BASE_COORD>1</BASE_COORD>
                    </READ_SPEC>
               </SPOT_DECODE_SPEC>
          </SPOT_DESCRIPTOR>
     </DESIGN>
     <PLATFORM>
          <ILLUMINA>
               <INSTRUMENT_MODEL>Illumina HiSeq 2000</INSTRUMENT_MODEL>
          </ILLUMINA>
     </PLATFORM>
     <PROCESSING>
          <PIPELINE>
               <PIPE_SECTION>
                    <STEP_INDEX>1</STEP_INDEX>
                    <PREV_STEP_INDEX>NIL</PREV_STEP_INDEX>
                    <PROGRAM/>
                    <VERSION/>
               </PIPE_SECTION>
          </PIPELINE>
     </PROCESSING>
     <EXPERIMENT_LINKS>
          <EXPERIMENT_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>DRS001521</ID>
               </XREF_LINK>
          </EXPERIMENT_LINK>
          <EXPERIMENT_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>DRR002116</ID>
               </XREF_LINK>
          </EXPERIMENT_LINK>
          <EXPERIMENT_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>DRA000567</ID>
               </XREF_LINK>
          </EXPERIMENT_LINK>
          <EXPERIMENT_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRX001563&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </EXPERIMENT_LINK>
          <EXPERIMENT_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRX001563&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </EXPERIMENT_LINK>
     </EXPERIMENT_LINKS>
     <EXPERIMENT_ATTRIBUTES>
          <EXPERIMENT_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>32568360</VALUE>
          </EXPERIMENT_ATTRIBUTE>
          <EXPERIMENT_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>3256836000</VALUE>
          </EXPERIMENT_ATTRIBUTE>
     </EXPERIMENT_ATTRIBUTES>
</EXPERIMENT>
</ROOT>"""


RUN_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="DRR002116&amp;display=xml">
<RUN run_date="2011-09-01T00:00:00+09:00" run_center="RIKEN_CDB" alias="DRR002116" center_name="RIKEN_CDB" accession="DRR002116" broker_name="DDBJ">
     <IDENTIFIERS>
          <PRIMARY_ID>DRR002116</PRIMARY_ID>
          <SUBMITTER_ID namespace="RIKEN_CDB">DRR002116</SUBMITTER_ID>
     </IDENTIFIERS>
     <TITLE>Illumina HiSeq 2000 sequencing; Exp_Gg_HH16_1_embryo_mRNAseq</TITLE>
     <EXPERIMENT_REF refname="DRX001563" refcenter="RIKEN_CDB" accession="DRX001563"/>
     <RUN_LINKS>
          <RUN_LINK>
               <XREF_LINK>
                    <DB>ENA-STUDY</DB>
                    <ID>DRP000595</ID>
               </XREF_LINK>
          </RUN_LINK>
          <RUN_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>DRS001521</ID>
               </XREF_LINK>
          </RUN_LINK>
          <RUN_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>DRA000567</ID>
               </XREF_LINK>
          </RUN_LINK>
          <RUN_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRR002116&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </RUN_LINK>
          <RUN_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRR002116&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </RUN_LINK>
     </RUN_LINKS>
     <RUN_ATTRIBUTES>
          <RUN_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>32568360</VALUE>
          </RUN_ATTRIBUTE>
          <RUN_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>3256836000</VALUE>
          </RUN_ATTRIBUTE>
          <RUN_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2013-07-19</VALUE>
          </RUN_ATTRIBUTE>
          <RUN_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2017-08-11</VALUE>
          </RUN_ATTRIBUTE>
     </RUN_ATTRIBUTES>
</RUN>
</ROOT>"""


SAMPLE_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="DRS001521&amp;display=xml">
<SAMPLE center_name="BioSample" alias="SAMD00004104" accession="DRS001521" broker_name="DDBJ">
     <IDENTIFIERS>
          <PRIMARY_ID>DRS001521</PRIMARY_ID>
          <EXTERNAL_ID namespace="BioSample">SAMD00004104</EXTERNAL_ID>
     </IDENTIFIERS>
     <TITLE>Gg_HH16_1_embryo_mRNAseq</TITLE>
     <SAMPLE_NAME>
          <TAXON_ID>9031</TAXON_ID>
          <SCIENTIFIC_NAME>Gallus gallus</SCIENTIFIC_NAME>
     </SAMPLE_NAME>
     <SAMPLE_LINKS>
          <SAMPLE_LINK>
               <XREF_LINK>
                    <DB>ENA-STUDY</DB>
                    <ID>DRP000595</ID>
               </XREF_LINK>
          </SAMPLE_LINK>
          <SAMPLE_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>DRX001563</ID>
               </XREF_LINK>
          </SAMPLE_LINK>
          <SAMPLE_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>DRR002116</ID>
               </XREF_LINK>
          </SAMPLE_LINK>
          <SAMPLE_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>DRA000567</ID>
               </XREF_LINK>
          </SAMPLE_LINK>
          <SAMPLE_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRS001521&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </SAMPLE_LINK>
          <SAMPLE_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRS001521&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </SAMPLE_LINK>
     </SAMPLE_LINKS>
     <SAMPLE_ATTRIBUTES>
          <SAMPLE_ATTRIBUTE>
               <TAG>sample_name</TAG>
               <VALUE>DRS001521</VALUE>
          </SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE>
               <TAG>sample comment</TAG>
               <VALUE>mRNAseq of chicken at stage HH16 (biological replicate 1)</VALUE>
          </SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>32568360</VALUE>
          </SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>3256836000</VALUE>
          </SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2013-07-20</VALUE>
          </SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2015-08-24</VALUE>
          </SAMPLE_ATTRIBUTE>
     </SAMPLE_ATTRIBUTES>
</SAMPLE>
</ROOT>"""


STUDY_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="DRP000595&amp;display=xml">
<STUDY accession="DRP000595" center_name="RIKEN_CDB" alias="DRP000595" broker_name="DDBJ">
     <IDENTIFIERS>
          <PRIMARY_ID>DRP000595</PRIMARY_ID>
          <SECONDARY_ID>PRJDB90</SECONDARY_ID>
          <EXTERNAL_ID namespace="BioProject" label="BioProject ID">PRJDB90</EXTERNAL_ID>
          <SUBMITTER_ID namespace="RIKEN_CDB">DRP000595</SUBMITTER_ID>
     </IDENTIFIERS>
     <DESCRIPTOR>
          <STUDY_TITLE>Whole transcriptome identification of turtle</STUDY_TITLE>
          <STUDY_TYPE existing_study_type="Transcriptome Analysis"/>
          <STUDY_ABSTRACT>Whole transcriptome of turtle (Pelodiscus sinensis) was identified in this study by using stranded sequencing methods.</STUDY_ABSTRACT>
          <CENTER_PROJECT_NAME>RIKEN_CDB</CENTER_PROJECT_NAME>
     </DESCRIPTOR>
     <STUDY_LINKS>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>pubmed</DB>
                    <ID>23624526</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>DRS001495-DRS001536</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>DRX001537-DRX001578</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>DRR002090-DRR002131</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMISSION</DB>
                    <ID>DRA000567</ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRP000595&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
          <STUDY_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRP000595&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </STUDY_LINK>
     </STUDY_LINKS>
     <STUDY_ATTRIBUTES>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>1371813555</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>158881910957</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-FIRST-PUBLIC</TAG>
               <VALUE>2013-07-19</VALUE>
          </STUDY_ATTRIBUTE>
          <STUDY_ATTRIBUTE>
               <TAG>ENA-LAST-UPDATE</TAG>
               <VALUE>2015-06-22</VALUE>
          </STUDY_ATTRIBUTE>
     </STUDY_ATTRIBUTES>
</STUDY>
</ROOT>"""


SUBMISSION_XML = """<?xml version="1.0" encoding="UTF-8"?>
<ROOT request="DRA000567&amp;display=xml">
<SUBMISSION accession="DRA000567" center_name="RIKEN_CDB" alias="DRA000567" lab_name="Group for Morphological Evolution, Center for Developmental Biology, Kobe Institute, RIKEN" submission_comment="Time course gene expression profiles of turtle (Pelodiscus sinensis) and chicken (Gallus gallus) embryos were examined. Whole transcriptome of turtle was also determined by uding stranded sequencing methods." broker_name="DDBJ">
     <IDENTIFIERS>
          <PRIMARY_ID>DRA000567</PRIMARY_ID>
          <SUBMITTER_ID namespace="RIKEN_CDB">DRA000567</SUBMITTER_ID>
     </IDENTIFIERS>
     <TITLE>Submitted by RIKEN_CDB on 19-JUL-2013</TITLE>
     <SUBMISSION_LINKS>
          <SUBMISSION_LINK>
               <XREF_LINK>
                    <DB>ENA-STUDY</DB>
                    <ID>DRP000593-DRP000595</ID>
               </XREF_LINK>
          </SUBMISSION_LINK>
          <SUBMISSION_LINK>
               <XREF_LINK>
                    <DB>ENA-SAMPLE</DB>
                    <ID>DRS001495-DRS001536</ID>
               </XREF_LINK>
          </SUBMISSION_LINK>
          <SUBMISSION_LINK>
               <XREF_LINK>
                    <DB>ENA-EXPERIMENT</DB>
                    <ID>DRX001537-DRX001578</ID>
               </XREF_LINK>
          </SUBMISSION_LINK>
          <SUBMISSION_LINK>
               <XREF_LINK>
                    <DB>ENA-RUN</DB>
                    <ID>DRR002090-DRR002131</ID>
               </XREF_LINK>
          </SUBMISSION_LINK>
          <SUBMISSION_LINK>
               <XREF_LINK>
                    <DB>ENA-FASTQ-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRA000567&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes]]></ID>
               </XREF_LINK>
          </SUBMISSION_LINK>
          <SUBMISSION_LINK>
               <XREF_LINK>
                    <DB>ENA-SUBMITTED-FILES</DB>
                    <ID><![CDATA[http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=DRA000567&result=read_run&fields=run_accession,submitted_ftp,submitted_md5,submitted_bytes,submitted_format]]></ID>
               </XREF_LINK>
          </SUBMISSION_LINK>
     </SUBMISSION_LINKS>
     <SUBMISSION_ATTRIBUTES>
          <SUBMISSION_ATTRIBUTE>
               <TAG>ENA-SPOT-COUNT</TAG>
               <VALUE>1371813555</VALUE>
          </SUBMISSION_ATTRIBUTE>
          <SUBMISSION_ATTRIBUTE>
               <TAG>ENA-BASE-COUNT</TAG>
               <VALUE>158881910957</VALUE>
          </SUBMISSION_ATTRIBUTE>
     </SUBMISSION_ATTRIBUTES>
</SUBMISSION>
</ROOT>"""
