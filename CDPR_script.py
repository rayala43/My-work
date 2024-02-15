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
    vcf_file = r'C:/Users/GenepoweRx_Madhu/Downloads/vcf_files_all/KHCDPRGPTTL126_final.vcf'
    output_file = r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/KHCDPRGPTTL126_final.vcf'

    bed_positions = read_bed_file(bed_file)
    filtered_vcf_records = filter_vcf_file(vcf_file, bed_positions)
    write_filtered_vcf(filtered_vcf_records, output_file)

if __name__ == "__main__":
    main()


vcf = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/KHCDPRGPTTL126_final.vcf', comment= '#', sep = '\t', header=None, low_memory=False)
vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

sample_cols = vcf['SAMPLE'].str.split(':', expand=True)
sample_cols.columns = ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']

# Assign the values to the newly created columns
vcf = pd.concat([vcf, sample_cols], axis=1)
vcf = vcf[['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL','RDF', 'RDR', 'ADF', 'ADR']]


vcf['HET'] = vcf['INFO'].str.extract(r'HET=(\d)')
vcf['HOM'] = vcf['INFO'].str.extract(r'HOM=(\d)')

# Create a new column 'Zygosity' based on conditions
vcf['Zygosity'] = ''

vcf.loc[vcf['HOM'] == '1', 'Zygosity'] = 'Homozygous'
vcf.loc[vcf['HET'] == '1', 'Zygosity'] = 'Heterozygous'
vcf['GT'] = vcf['GT'].astype(str)

vcf["Gene Name"] = vcf["INFO"].str.extract('GENEINFO=(?P<GENEINFO>.+?);')
vcf['Gene Name'] = vcf['Gene Name'].apply(lambda x: ','.join(set([segment.split(':')[0] for segment in x.split('|')])) if pd.notnull(x) else '')

vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
vcf['csq'] = vcf['CSQ'].str.split(',')
vcf = vcf.explode('csq')

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

vcf['Protein Position and Amino Acid'] = vcf['Amino_acids'].str[0] + vcf['Protein_position'] + np.where(vcf['Amino_acids'].str[-1] == vcf['Amino_acids'].str[0], '', vcf['Amino_acids'].str[-1])

vcf[['HGVSc', 'HGVSc (Transcript)']] = vcf['HGVSC'].str.split(':', 1, expand=True)
vcf[['HGVSp', 'HGVSp (Transcript)']] = vcf['HGVSP'].str.split(':', 1, expand=True)

vcf_final = vcf[['Gene Name', 'rsID','CHROM', 'POS', 'REF', 'ALT', 'Zygosity', 'Consequence', 'IMPACT',
          'ClinVar_CLNDN', 'CLIN_SIG', 'ClinVar_CLNREVSTAT',
          'ClinVar', 'HGVSc', 'HGVSc (Transcript)', 'HGVSp', 'HGVSp (Transcript)', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RDF', 'RDR', 'ADF',
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

vcf_final['consequence'] = vcf_final['Consequence'].str.split(',').str[0]

df_1 = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Madhu_folder_04_07_2023/kidney_health_final.vcf/consequence.xlsx')

merged_1 = pd.merge(vcf_final, df_1, on='consequence', how='left', sort=False)

df_2 = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Madhu_folder_04_07_2023/kidney_health_final.vcf/IMPACT.xlsx')

merged_2 = pd.merge(merged_1, df_2, on = 'IMPACT', how='left', sort=False)

merged_2 = merged_2[['Gene Name', 'rsID', 'CHROM', 'POS', 'REF', 'ALT', 'Zygosity',
       'Consequence','Consequence_score', 'IMPACT', 'IMPACT_score', 'ClinVar_CLNDN', 'CLIN_SIG',
       'ClinVar_CLNREVSTAT', 'ClinVar', 'HGVSc', 'HGVSc (Transcript)', 'HGVSp',
       'HGVSp (Transcript)', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ',
       'PVAL', 'RDF', 'RDR', 'ADF', 'ADR', 'SIFT', 'PolyPhen', 'AF', 'AFR_AF',
       'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomADe_AF', 'gnomADe_AFR_AF',
       'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF',
       'gnomADe_NFE_AF', 'gnomADe_OTH_AF', 'gnomADe_SAS_AF', 'gnomADg_AF',
       'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF',
       'gnomADg_EAS_AF', 'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF',
       'gnomADg_OTH_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'BIOTYPE',
       'EXON', 'INTRON', 'Protein Position and Amino Acid', 'Codons', 'STRAND',
       'PUBMED']]

df_gene = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/OneDrive/Desktop/cardiac_genes.xlsx')

merged_2['Gene_Match'] = 'No'

# Iterate through each gene in df1
for genes in merged_2['Gene Name']:
    if isinstance(genes, str):  # Check if the gene value is a non-null string
        gene_list = genes.split(',')  # Split the genes by comma to create a list
        match = any(gene in df_gene['Gene Name'].values for gene in gene_list)  # Check if any gene in the list exists in df2
        if match:
            merged_2.loc[merged_2['Gene Name'] == genes, 'Gene_Match'] = 'Yes'
            
merged_2 = merged_2[['Gene Name', 'Gene_Match', 'rsID', 'CHROM', 'POS', 'REF', 'ALT', 'Zygosity',
       'Consequence', 'Consequence_score', 'IMPACT', 'IMPACT_score',
       'ClinVar_CLNDN', 'CLIN_SIG', 'ClinVar_CLNREVSTAT', 'ClinVar', 'HGVSc',
       'HGVSc (Transcript)', 'HGVSp', 'HGVSp (Transcript)', 'GT', 'GQ', 'SDP',
       'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RDF', 'RDR', 'ADF', 'ADR', 'SIFT',
       'PolyPhen', 'AF', 'AFR_AF',
       'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomADe_AF', 'gnomADe_AFR_AF',
       'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF',
       'gnomADe_NFE_AF', 'gnomADe_OTH_AF', 'gnomADe_SAS_AF', 'gnomADg_AF',
       'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF',
       'gnomADg_EAS_AF', 'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF',
       'gnomADg_OTH_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'BIOTYPE', 'EXON', 'INTRON',
       'Protein Position and Amino Acid', 'Codons', 'STRAND', 'PUBMED']]
merged_2['POS'] = merged_2['POS'].astype('int64')

import pandas as pd
x = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/OneDrive/Desktop/Cardiac_Lit_final_hg38_hg37.xlsx')
x['chrom'] = x['Chrom-pos-Ref-Alt_38'].str.split(',')
x = x.explode('chrom')

x['CHROM'] = x['chrom'].str.split('-').str[0]

# Function to add 'chr' prefix conditionally
def add_chr_prefix(chrom):
    if pd.notnull(chrom) and chrom.strip() != '':
        return 'chr' + str(chrom)
    else:
        return chrom

# Applying the function to the 'chromosome' column
x['CHROM'] = x['CHROM'].apply(add_chr_prefix)
x['CHROM'] = x['CHROM'].str.strip()
x['CHROM'] = x['CHROM'].str.replace(r'\s+', '')
x['POS'] = x['chrom'].str.split('-').str[1]
x['REF'] = x['chrom'].str.split('-').str[2]
x['ALT'] = x['chrom'].str.split('-').str[3]


x.dropna(subset=['CHROM'], inplace=True)
# Drop rows with empty cells after removing leading and trailing whitespaces
x['CHROM'] = x['CHROM'].str.strip()
x['POS'] = x['POS'].str.strip()
# Dropping rows with empty cells and NaN values in both 'chromosome' and 'position' columns
x.dropna(subset=['CHROM', 'POS'], inplace=True)
df_3 = x[['CHROM', 'POS', 'REF', 'ALT']]
df_3['Literature'] = 'Yes'
df_3.drop_duplicates(subset='POS', inplace=True)
df_3['POS'] = df_3['POS'].astype('int64')
df_3 = df_3.reset_index()
df_3 = df_3[['CHROM', 'POS', 'REF', 'ALT', 'Literature']]

df = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/KAPA HyperExome_hg38_capture_targets (1).bed', sep = '\t', header = None)
df.columns = ['chromosome', 'Start_pos', 'End_pos', 'INFO']

df['Extended_Start_pos'] = df['Start_pos'] - 20
df['Extended_End_pos'] = df['End_pos'] + 20

df['gene_symbol'] = df['INFO'].str.extract(r'gene_symbol=([^;]+)')
df = df[['chromosome', 'Extended_Start_pos', 'Extended_End_pos', 'INFO', 'gene_symbol']]


# Step 1: Create a dictionary from the df DataFrame
chromosome_dict = {}
for _, row in df.iterrows():
    chromosome = row['chromosome']
    start_pos = row['Extended_Start_pos']
    end_pos = row['Extended_End_pos']
    if chromosome not in chromosome_dict:
        chromosome_dict[chromosome] = []
    chromosome_dict[chromosome].append((start_pos, end_pos))

# Step 2: Define a function to check coverage
def check_coverage(row):
    pos = row['POS']
    chromosome = row['CHROM']
    if chromosome in chromosome_dict:
        ranges = chromosome_dict[chromosome]
        for start, end in ranges:
            if start <= pos <= end:
                return 'Covered'
    return 'Not_Covered'

# Step 3: Apply the function to create the new column in dataset2
df_3['Covered/Not_Covered'] = df_3.apply(check_coverage, axis=1)

df_3 = df_3[df_3['Covered/Not_Covered'] == 'Covered']

merged_2['POS'] = merged_2['POS'].astype('int64')
df_3['POS'] = df_3['POS'].astype('int64')
merged_3 = pd.merge(merged_2, df_3, on=['CHROM', 'POS', 'REF', 'ALT'], how='left', sort=False)
merged_3['Literature'] = merged_3['Literature'].fillna('No')

merged_3 = merged_3[['Gene Name', 'Gene_Match', 'rsID', 'Literature', 'CHROM', 'POS', 'REF', 'ALT', 'Zygosity',
       'Consequence','Consequence_score', 'IMPACT', 'IMPACT_score', 'ClinVar_CLNDN', 'CLIN_SIG',
       'ClinVar_CLNREVSTAT', 'ClinVar', 'HGVSc', 'HGVSc (Transcript)', 'HGVSp',
       'HGVSp (Transcript)', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ',
       'PVAL', 'RDF', 'RDR', 'ADF', 'ADR', 'SIFT', 'PolyPhen', 'AF', 'AFR_AF',
       'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'gnomADe_AF', 'gnomADe_AFR_AF',
       'gnomADe_AMR_AF', 'gnomADe_ASJ_AF', 'gnomADe_EAS_AF', 'gnomADe_FIN_AF',
       'gnomADe_NFE_AF', 'gnomADe_OTH_AF', 'gnomADe_SAS_AF', 'gnomADg_AF',
       'gnomADg_AFR_AF', 'gnomADg_AMI_AF', 'gnomADg_AMR_AF', 'gnomADg_ASJ_AF',
       'gnomADg_EAS_AF', 'gnomADg_FIN_AF', 'gnomADg_MID_AF', 'gnomADg_NFE_AF',
       'gnomADg_OTH_AF', 'gnomADg_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'BIOTYPE',
       'EXON', 'INTRON', 'Protein Position and Amino Acid', 'Codons', 'STRAND',
       'PUBMED']]

merged_3.to_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Processed_vcf_files/KHCDPRGPTTL126_final_depth_vcf_processed.xlsx', index=False)