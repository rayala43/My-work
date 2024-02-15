df = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/1232_PGx_unique_star_alleles_Function_Evidence.xlsx')
df = df[['CHROM', 'POS', 'REF', 'ALT', 'rsID', 'Gene', 'Gene_star_allele_updated']]

data = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/KHAIGPRX1_snp_pos_distance.xlsx')
data['HET'] = data['INFO'].str.extract(r'HET=(\d)')
data['HOM'] = data['INFO'].str.extract(r'HOM=(\d)')
# Create a new column 'Zygosity' based on conditions
data['Zygosity'] = ''
data.loc[data['HOM'] == '1', 'Zygosity'] = 'Homozygous'
data.loc[data['HET'] == '1', 'Zygosity'] = 'Heterozygous'
data = data[data['Haplotype'] == 'Snp']
data = data[['CHROM', 'POS', 'REF', 'ALT', 'rsID', 'Covered/Not_Covered', 'Gene', 'Zygosity']]
# Function to find the nearest POS and corresponding 'Gene_star_allele_updated'
def find_nearest_pos(row):
    df2_subset = df[df['Gene'] == row['Gene']]
    if not df2_subset.empty:
        nearest_pos_index = (df2_subset['POS'] - row['POS']).abs().idxmin()
        nearest_pos = df2_subset.loc[nearest_pos_index, 'POS']
        nearest_gene_star_allele_updated = df2_subset.loc[nearest_pos_index, 'Gene_star_allele_updated']
        return nearest_pos, nearest_gene_star_allele_updated
    else:
        return None, None

# Apply the function to create 'POS_df2' and 'Gene_star_allele_df2' columns in data
data[['POS_df2', 'Gene_star_allele_df2']] = data.apply(find_nearest_pos, axis=1, result_type='expand')
data['POS_diff'] = data['POS_df2'] - data['POS']
data = data.dropna(subset=['POS_df2'])
data