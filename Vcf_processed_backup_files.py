import numpy as np
import pandas as pd
import polars as pl
import sys
import re
import os
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import plotly.express as px
import warnings
warnings.filterwarnings("ignore")
pd.set_option('display.max_columns',None)
import psycopg2

print('Libraries Loaded')

def read_bed_file(bed_file):
    bed_positions = set()
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines if present
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom = fields[0]
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    continue  # Skip this line if start or end position is not an integer
                for pos in range(start, end + 1):
                    bed_positions.add((chrom, pos))
    return bed_positions

def normalize_chrom_name(chrom):
    return chrom.split('_')[0]

def filter_vcf_file(vcf_file, bed_positions):
    filtered_vcf_records = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Preserve header lines in the output
                filtered_vcf_records.append(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                raw_chrom = fields[0]
                chrom = normalize_chrom_name(raw_chrom)
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue  # Skip this line if 'POS' is not an integer
                if (chrom, pos) in bed_positions:
                    filtered_vcf_records.append(line)
    return filtered_vcf_records

def write_filtered_vcf(filtered_vcf_records, output_file):
    with open(output_file, 'w') as f:
        for record in filtered_vcf_records:
            f.write(record)

def main():
    bed_file = r'C:/Users/GenepoweRx_Madhu/Downloads/BED_files/srinivas_sir_covered.bed'
    vcf_file = r'C:/Users/GenepoweRx_Madhu/Downloads/vcf_files_all/KHCDPRGPTTL124_final.vcf'
    output_file = r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/KHCDPRGPTTL124_final.vcf'

    bed_positions = read_bed_file(bed_file)
    filtered_vcf_records = filter_vcf_file(vcf_file, bed_positions)
    write_filtered_vcf(filtered_vcf_records, output_file)

if __name__ == "__main__":
    main()

print('Covered variants Extracted from the vcf file without Distribung the Metadata')

vcf = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/KHCDPRGPTTL124_final.vcf', comment= '#', sep = '\t', header=None, low_memory=False)
vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

sample_cols = vcf['SAMPLE'].str.split(':', expand=True)
sample_cols.columns = ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'Rvcf', 'RDR', 'Avcf', 'ADR']

# Assign the values to the newly created columns
vcf = pd.concat([vcf, sample_cols], axis=1)
vcf = vcf[['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL','Rvcf', 'RDR', 'Avcf', 'ADR']]

print('Covered Variants Data Loaded')

vcf['HET'] = vcf['INFO'].str.extract(r'HET=(\d)')
vcf['HOM'] = vcf['INFO'].str.extract(r'HOM=(\d)')

# Create a new column 'Zygosity' based on conditions
vcf['Zygosity'] = ''

vcf.loc[vcf['HOM'] == '1', 'Zygosity'] = 'Homozygous'
vcf.loc[vcf['HET'] == '1', 'Zygosity'] = 'Heterozygous'
vcf['GT'] = vcf['GT'].astype(str)

print('Zygosity Column Extracted')

vcf["Gene Name"] = vcf["INFO"].str.extract('GENEINFO=(?P<GENEINFO>.+?);')
vcf['Gene Name'] = vcf['Gene Name'].apply(lambda x: ','.join(set([segment.split(':')[0] for segment in x.split('|')])) if pd.notnull(x) else '')

vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
vcf['csq'] = vcf['CSQ'].str.split(',')
vcf = vcf.explode('csq')

print('Genes and CSQ columns splitted')

########################################################### Required columns extraction from the CSQ column ####################
vcf['ClinVar_CLNDN'] = vcf['csq'].str.split('|').str[82]
vcf['CLIN_SIG'] = vcf['csq'].str.split('|').str[70]
vcf['ClinVar_CLNREVSTAT'] = vcf['csq'].str.split('|').str[81]
vcf['ClinVar'] = vcf['csq'].str.split('|').str[79]
vcf['HGVSC'] = vcf['csq'].str.split('|').str[10]
vcf['HGVSP'] = vcf['csq'].str.split('|').str[11]
vcf['PolyPhen'] = vcf['csq'].str.split('|').str[38]
vcf['BIOTYPE'] = vcf['csq'].str.split('|').str[7]
vcf['EXON'] = vcf['csq'].str.split('|').str[8]
vcf['INTRON'] = vcf['csq'].str.split('|').str[9]
vcf['Protein_position'] = vcf['csq'].str.split('|').str[14]
vcf['Amino_acids'] = vcf['csq'].str.split('|').str[15]
vcf['Codons'] = vcf['csq'].str.split('|').str[16]
vcf['STRAND'] = vcf['csq'].str.split('|').str[19]
vcf['PUBMED'] = vcf['csq'].str.split('|').str[73]
vcf['Consequence'] = vcf['csq'].str.split('|').str[1]
vcf['IMPACT'] = vcf['csq'].str.split('|').str[2]
vcf['SIFT'] = vcf['csq'].str.split('|').str[37]
################################################## Frequency columns extraction ################################################
vcf['AF'] = vcf['csq'].str.split('|').str[42]
vcf['AFR_AF'] = vcf['csq'].str.split('|').str[43]
vcf['AMR_AF'] = vcf['csq'].str.split('|').str[44]
vcf['EAS_AF'] = vcf['csq'].str.split('|').str[45]
vcf['EUR_AF'] = vcf['csq'].str.split('|').str[46]
vcf['SAS_AF'] = vcf['csq'].str.split('|').str[47]
vcf['gnomADe_AF'] = vcf['csq'].str.split('|').str[48]
vcf['gnomADe_AFR_AF'] = vcf['csq'].str.split('|').str[49]
vcf['gnomADe_AMR_AF'] = vcf['csq'].str.split('|').str[50]
vcf['gnomADe_ASJ_AF'] = vcf['csq'].str.split('|').str[51]
vcf['gnomADe_EAS_AF'] = vcf['csq'].str.split('|').str[52]
vcf['gnomADe_FIN_AF'] = vcf['csq'].str.split('|').str[53]
vcf['gnomADe_NFE_AF'] = vcf['csq'].str.split('|').str[54]
vcf['gnomADe_OTH_AF'] = vcf['csq'].str.split('|').str[55]
vcf['gnomADe_SAS_AF'] = vcf['csq'].str.split('|').str[56]
vcf['gnomADg_AF'] = vcf['csq'].str.split('|').str[57]
vcf['gnomADg_AFR_AF'] = vcf['csq'].str.split('|').str[58]
vcf['gnomADg_AMI_AF'] = vcf['csq'].str.split('|').str[59]
vcf['gnomADg_AMR_AF'] = vcf['csq'].str.split('|').str[60]
vcf['gnomADg_ASJ_AF'] = vcf['csq'].str.split('|').str[61]
vcf['gnomADg_EAS_AF'] = vcf['csq'].str.split('|').str[62]
vcf['gnomADg_FIN_AF'] = vcf['csq'].str.split('|').str[63]
vcf['gnomADg_MID_AF'] = vcf['csq'].str.split('|').str[64]
vcf['gnomADg_NFE_AF'] = vcf['csq'].str.split('|').str[65]
vcf['gnomADg_OTH_AF'] = vcf['csq'].str.split('|').str[66]
vcf['gnomADg_SAS_AF'] = vcf['csq'].str.split('|').str[67]
vcf['MAX_AF'] = vcf['csq'].str.split('|').str[68]
vcf['MAX_AF_POPS'] = vcf['csq'].str.split('|').str[69]

print('Required Columns Extracted')

vcf['Protein Position and Amino Acid'] = vcf['Amino_acids'].str[0] + vcf['Protein_position'] + np.where(vcf['Amino_acids'].str[-1] == vcf['Amino_acids'].str[0], '', vcf['Amino_acids'].str[-1])

vcf[['HGVSc', 'HGVSc (Transcript)']] = vcf['HGVSC'].str.split(':', 1, expand=True)
vcf[['HGVSp', 'HGVSp (Transcript)']] = vcf['HGVSP'].str.split(':', 1, expand=True)

print('Protien and amino acid change and hgvsc, hgvsp transcripts extracted')

vcf_final = vcf[['Gene Name', 'rsID','CHROM', 'POS', 'REF', 'ALT', 'Zygosity', 'Consequence', 'IMPACT',
          'ClinVar_CLNDN', 'CLIN_SIG', 'ClinVar_CLNREVSTAT',
          'ClinVar', 'HGVSc', 'HGVSc (Transcript)', 'HGVSp', 'HGVSp (Transcript)', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'Rvcf', 'RDR', 'Avcf',
       'ADR', 'SIFT', 'PolyPhen', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF',
       'gnomADe_AF', 'gnomADe_AFR_AF', 'gnomADe_AMR_AF', 'gnomADe_ASJ_AF',
       'gnomADe_EAS_AF', 'gnomADe_FIN_AF', 'gnomADe_NFE_AF', 'gnomADe_OTH_AF',
       'gnomADe_SAS_AF', 'gnomADg_AF', 'gnomADg_AFR_AF', 'gnomADg_AMI_AF',
       'gnomADg_AMR_AF', 'gnomADg_ASJ_AF', 'gnomADg_EAS_AF', 'gnomADg_FIN_AF',
       'gnomADg_MID_AF', 'gnomADg_NFE_AF', 'gnomADg_OTH_AF', 'gnomADg_SAS_AF',
       'MAX_AF', 'MAX_AF_POPS', 'BIOTYPE', 'EXON', 'INTRON',
          'Protein Position and Amino Acid', 'Codons', 'STRAND', 'PUBMED']]

# Define the terms to remove
remove_terms = set(["not_specified", "not_provided"])

# Apply the filtering operation to 'Column1' only
vcf_final['ClinVar_CLNDN'] = vcf_final['ClinVar_CLNDN'].apply(lambda row: "&".join(
    [term for term in row.split("&") if term not in remove_terms]
    ) if isinstance(row, str) and not all(term in remove_terms for term in row.split("&")) else row)


vcf_final['CLIN_SIG'] = vcf_final['CLIN_SIG'].apply(lambda row: "&".join(
    [term for term in row.split("&") if term not in remove_terms]
    ) if isinstance(row, str) and not all(term in remove_terms for term in row.split("&")) else row)


vcf_final['ClinVar_CLNREVSTAT'] = vcf_final['ClinVar_CLNREVSTAT'].apply(lambda row: "&".join(
    [term for term in row.split("&") if term not in remove_terms]
    ) if isinstance(row, str) and not all(term in remove_terms for term in row.split("&")) else row)

vcf_final = vcf_final.astype(str).applymap(lambda x: x.replace('&', ',').replace('_', ' '))

print('Removed the unnecessary strings and unnecessary characters from the columns')

vcf_final['consequence'] = vcf_final['Consequence'].str.split(',').str[0]

vcf_1 = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Madhu_folder_04_07_2023/kidney_health_final.vcf/consequence.xlsx')

merged_1 = pd.merge(vcf_final, vcf_1, on='consequence', how='left', sort=False)

vcf_2 = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Madhu_folder_04_07_2023/kidney_health_final.vcf/IMPACT.xlsx')

print('Consequence scores and Impact scores Mapped to the consequence and Impact columns')

merged_2 = pd.merge(merged_1, vcf_2, on = 'IMPACT', how='left', sort=False)

merged_2 = merged_2[['Gene Name', 'rsID', 'CHROM', 'POS', 'REF', 'ALT', 'Zygosity',
       'Consequence','Consequence_score', 'IMPACT', 'IMPACT_score', 'ClinVar_CLNDN', 'CLIN_SIG',
       'ClinVar_CLNREVSTAT', 'ClinVar', 'HGVSc', 'HGVSc (Transcript)', 'HGVSp',
       'HGVSp (Transcript)', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ',
       'PVAL', 'Rvcf', 'RDR', 'Avcf', 'ADR', 'SIFT', 'PolyPhen', 'AF', 'AFR_AF',
       'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomADe_AF', 'gnomADe_AFR_AF',
       'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF',
       'gnomADe_NFE_AF', 'gnomADe_OTH_AF', 'gnomADe_SAS_AF', 'gnomADg_AF',
       'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF',
       'gnomADg_EAS_AF', 'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF',
       'gnomADg_OTH_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'BIOTYPE',
       'EXON', 'INTRON', 'Protein Position and Amino Acid', 'Codons', 'STRAND',
       'PUBMED']]

print('Reordered the required columns and exporting the file as excel file')

merged_2.to_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Processed_vcf_files/Backup_files/KHCDPRGPTTL124_final.xlsx', index=False)

print('Final file Exported')