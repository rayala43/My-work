import numpy as np
import pandas as pd
pd.set_option('display.max_columns',None)
import psycopg2
import warnings
warnings.filterwarnings("ignore")

# Read the VCF file into a Pandas DataFrame (skip first 80 rows as shown in a previous answer)
vcf = pd.read_csv(
    r'C:/Users/GenepoweRx_Madhu/Downloads/718_final.vcf',
    comment='#',
    sep='\t',
    header=None,
    low_memory=False,
    skiprows=80
)

# Assign appropriate column names
vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
vcf["Gene"] = vcf["INFO"].str.extract('GENEINFO=(?P<GENEINFO>.+?);')
vcf['Gene'] = vcf['Gene'].apply(lambda x: ','.join(set([segment.split(':')[0] for segment in x.split('|')])) if pd.notnull(x) else '')
df_gene = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Desktop/selective_genes.xlsx')

vcf['Matched_gene'] = 'No'

# Iterate through each gene in df1
for genes in vcf['Gene']:
    if isinstance(genes, str):  # Check if the gene value is a non-null string
        gene_list = genes.split(',')  # Split the genes by comma to create a list
        match = any(gene in df_gene['Gene'].values for gene in gene_list)  # Check if any gene in the list exists in df2
        if match:
            vcf.loc[vcf['Gene'] == genes, 'Matched_gene'] = 'Yes'
            
vcf = vcf[vcf['Matched_gene'] == 'Yes']

sample_cols = vcf['SAMPLE'].str.split(':', expand=True)
sample_cols.columns = ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']

# Assign the values to the newly created columns
vcf = pd.concat([vcf, sample_cols], axis=1)

# Create empty columns
columns = ['ADP', 'WT', 'HET', 'HOM', 'NC', 'CDA', 'OTH', 'S3D', 'WTD', 'dbSNPBuildID', 'SLO',
           'NSF', 'R3', 'R5', 'NSN', 'NSM', 'G5A', 'COMMON', 'RS', 'RV', 'TPA', 'CFL', 'GNO',
           'VLD', 'ASP', 'ASS', 'Ref', 'U3', 'U5', 'TOPMED', 'WGT', 'MTP', 'LSD', 'NOC',
           'DSS', 'SYN', 'KGPhase3', 'CAF', 'VC', 'MUT', 'KGPhase1', 'NOV', 'VP', 'SAO',
           'GENEINFO', 'INT', 'G5', 'OM', 'PMC', 'SSR', 'RSPOS', 'HD', 'PM']

for col in columns:
    vcf[col] = ''

# Populate columns based on 'INFO' values
for i, row in vcf.iterrows():
    info = row['INFO']
    items = info.split(';')
    for item in items:
        key_value = item.split('=')
        key = key_value[0]
        if key in columns:
            if len(key_value) > 1:
                value = key_value[1]
                vcf.at[i, key] = f"{key}={value}"

# Extract the numeric part of the HET column and convert it to an integer
vcf['HET_new'] = vcf['HET'].str.extract('(\d+)').astype(int)

# Create the Zygosity column based on the HET column
vcf['Zygosity'] = vcf['HET_new'].apply(lambda x: 'Heterozygous' if x == 1 else 'Homozygous')

# Extract the CSQ field from the INFO column
vcf['CSQ'] = vcf['INFO'].str.extract('CSQ=([^;]+)')

# Split the CSQ field into individual transcripts
transcripts = vcf['CSQ'].str.split(',')

CSQ_columns = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "SOURCE", "ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN"]

# Initialize an empty DataFrame to store the extracted values
extracted_df = pd.DataFrame(columns=CSQ_columns)

# Helper function to remove leading and trailing commas
def remove_leading_trailing_commas(s):
    if isinstance(s, str):
        return s.strip(',')
    return s

# Iterate over each transcript and extract values into separate columns
for transcript in transcripts:
    values = {}
    for entry in transcript:
        entry = remove_leading_trailing_commas(entry)  # Remove leading and trailing commas
        parts = entry.split('|')
        for i, column_name in enumerate(CSQ_columns):
            if i < len(parts):
                if column_name in values:
                    values[column_name] += ',' + parts[i]
                else:
                    values[column_name] = parts[i]
    extracted_df = extracted_df.append(values, ignore_index=True)

# Apply set to each cell to get unique values
unique_extracted_df = extracted_df.applymap(lambda x: ','.join(set(x.split(',')) if isinstance(x, str) else ''))

# Remove leading and trailing commas in the entire DataFrame
unique_extracted_df = unique_extracted_df.applymap(remove_leading_trailing_commas)

# Reset the index of both DataFrames
vcf.reset_index(drop=True, inplace=True)
unique_extracted_df.reset_index(drop=True, inplace=True)
combined_df = pd.concat([vcf, unique_extracted_df], axis=1)
columns_to_drop = [7, 8, 9, 11, 79, 81]
combined_df = combined_df.drop(columns=combined_df.columns[columns_to_drop])
combined_df = combined_df.replace({'': 'NA', np.nan: 'NA'})
combined_df.to_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/biotype_final_csq_specific_gene.xlsx', index = False)
print(f"Genes present in VCF file have been saved to {combined_df}")


#CSQ_columns = [
#    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
#    "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
#    "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
#    "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL",
#    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT",
#    "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "SOURCE", "GENE_PHENO", "SIFT", "PolyPhen",
#    "DOMAINS", "miRNA", "HGVS_OFFSET", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF",
#    "SAS_AF", "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF",
#    "gnomADe_EAS_AF", "gnomADe_FIN_AF", "gnomADe_NFE_AF", "gnomADe_OTH_AF", "gnomADe_SAS_AF",
#    "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF",
#    "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_OTH_AF",
#    "gnomADg_SAS_AF", "MAX_AF", "MAX_AF_POPS", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED",
#    "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS",
#    "ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN"]