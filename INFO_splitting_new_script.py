import numpy as np
import pandas as pd
import polars as pl
import sys
import re
import os
pd.set_option('display.max_columns',None)
import warnings
warnings.filterwarnings("ignore")

vcf = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/KHGLBS598_vep.vcf', comment= '#', sep = '\t', header=None, low_memory=False)
vcf.columns = ['CHROM', 'POS', 'rsID', 'Ref', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
sample_cols = vcf['SAMPLE'].str.split(':', expand=True)
sample_cols.columns = ['GT', 'SO', 'REPCN', 'REPCI', 'ADSP', 'ADFL', 'ADIR', 'LC']

# Assign the values to the newly created columns
vcf = pd.concat([vcf, sample_cols], axis=1)
vcf = vcf[['CHROM', 'POS', 'rsID', 'Ref', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT', 'SO', 'REPCN', 'REPCI', 'ADSP', 'ADFL', 'ADIR', 'LC']]
vcf_1 = vcf.copy()

# Create empty columns
columns = ['END', 'REF', 'REPID', 'RL', 'RU', 'SVTYPE', 'VARID']

for col in columns:
    vcf_1[col] = ''

# Populate columns based on 'info' values
for i, row in vcf_1.iterrows():
    info = row['INFO']
    items = info.split(';')
    for item in items:
        key_value = item.split('=')
        key = key_value[0]
        if key in columns:
            if len(key_value) > 1:
                value = key_value[1]
                vcf_1.at[i, key] = value  # Extract the value after '='
            else:
                vcf_1.at[i, key] = None  # Handle the case where there's no value
        else:
            vcf_1.at[i, key] = 'null'

column_names=['Allele','Consequence','IMPACT','SYMBOL','Gene','Feature_type','Feature','BIOTYPE','EXON','INTRON','HGVSc','HGVSp','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons','Existing_variation','DISTANCE','STRAND','FLAGS','VARIANT_CLASS','SYMBOL_SOURCE','HGNC_ID','CANONICAL','MANE_SELECT','MANE_PLUS_CLINICAL','TSL','APPRIS','CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','UNIPROT_ISOFORM','SOURCE','GENE_PHENO','SIFT','PolyPhen','DOMAINS','miRNA','HGVS_OFFSET','AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF','gnomADe_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_ASJ_AF','gnomADe_EAS_AF','gnomADe_FIN_AF','gnomADe_NFE_AF','gnomADe_OTH_AF','gnomADe_SAS_AF','gnomADg_AF','gnomADg_AFR_AF','gnomADg_AMI_AF','gnomADg_AMR_AF','gnomADg_ASJ_AF','gnomADg_EAS_AF','gnomADg_FIN_AF','gnomADg_MID_AF','gnomADg_NFE_AF','gnomADg_OTH_AF','gnomADg_SAS_AF','MAX_AF','MAX_AF_POPS','CLIN_SIG','SOMATIC','PHENO','PUBMED','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','TRANSCRIPTION_FACTORS','ClinVar','ClinVar_CLNSIG','ClinVar_CLNREVSTAT','ClinVar_CLNDN']

def process_csq(column_index, csq):
    values = [transcript.split('|')[column_index] for transcript in csq if transcript.split('|')[column_index]]
    return ','.join(set(values))

# Iterate through column names and apply processing function to create new columns
for i, col in enumerate(column_names):
    vcf_1[col] = vcf_1['INFO'].str.extract(r'CSQ=(.*)')
    vcf_1[col] = vcf_1[col].str.split(',').apply(lambda csq: process_csq(i, csq))
print("CSQ_splitting is done: ")
vcf_1 = vcf_1.drop(['INFO', 'CSQ'], axis=1)
vcf_1.to_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/KHGLBS598_vep_splitted.xlsx', index = False)
vcf_1