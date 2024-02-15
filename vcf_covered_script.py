import pandas as pd

vcf = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/vcf_files_all/KHAIGPRX1351_final.vcf', comment= '#', sep = '\t', header=None, low_memory=False)
vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

sample_cols = vcf['SAMPLE'].str.split(':', expand=True)
sample_cols.columns = ['GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RBQ', 'ABQ', 'RDF', 'RDR', 'ADF', 'ADR']

# Assign the values to the newly created columns
vcf = pd.concat([vcf, sample_cols], axis=1)
vcf = vcf[['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT', 'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL','RDF', 'RDR', 'ADF', 'ADR']]
vcf['DP'] = vcf['DP'].astype('int')
vcf = vcf[vcf['DP'] >= 15]

df = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/KAPA HyperExome_hg38_capture_targets (1).bed', sep = '\t', header = None, error_bad_lines=False)
df.columns = ['chromosome', 'Start_pos', 'End_pos', 'INFO']
df['Extended_Start_pos'] = df['Start_pos'] - 20
df['Extended_End_pos'] = df['End_pos'] + 20
df = df[['chromosome', 'Extended_Start_pos', 'Extended_End_pos']]

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

# Step 3: Apply the function to create the new column in data
vcf['Covered/Not_Covered'] = vcf.apply(check_coverage, axis=1)
#vcf = vcf[vcf['Covered/Not_Covered'] == 'Covered']
vcf["Gene Name"] = vcf["INFO"].str.extract('GENEINFO=(?P<GENEINFO>.+?);')
vcf['Gene Name'] = vcf['Gene Name'].apply(lambda x: ','.join(set([segment.split(':')[0] for segment in x.split('|')])) if pd.notnull(x) else '')
#vcf.to_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/KHAIGPRX1331_snp_pos.xlsx', index=False)
vcf