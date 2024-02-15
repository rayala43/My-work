import pandas as pd
import re
import sys
df = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/KHAIGPRX959_final_DP.vcf', comment= '#', sep = '\t', header=None, low_memory=False)
df.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
df['HET'] = df['INFO'].str.extract(r'HET=(\d)')
df['HOM'] = df['INFO'].str.extract(r'HOM=(\d)')
# Create a new column 'Zygosity' based on conditions
df['Zygosity'] = ''
df.loc[df['HOM'] == '1', 'Zygosity'] = 'Homozygous'
df.loc[df['HET'] == '1', 'Zygosity'] = 'Heterozygous'
df = df.drop(columns=['HET', 'HOM'], axis = 1)
data = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/RSIDS_POSITION_DBSNP.xlsx')
df = pd.merge(df, data, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'inner', sort = False)
df2 = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/13_Genes_Key_Variants_Coordinates.csv', sep = '\t')
df2 = df2[['CHROM', 'POS', 'REF', 'ALT', 'Haplotype']]
merged = pd.merge(df, df2, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'left', sort = False)
madhu = merged.copy()
madhu['rsID_mapped'] = madhu.groupby('Haplotype')['rsID'].transform(lambda x: ','.join(set(x)))
madhu['Zygosity_mapped'] = madhu.groupby(['Haplotype'])['Zygosity'].transform(lambda x: ','.join(x))
madhu.drop_duplicates(subset='Haplotype', inplace=True)
madhu['rsID_mapped_sorted'] = madhu['rsID_mapped'].apply(lambda x: ','.join(sorted(x.split(','))))
madhu = madhu.drop(columns=['rsID_mapped', 'rsID_updated', 'Haplotype'], axis=1)
data = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Multiple_Positions.xlsx')
data = data.rename({'Haplotype':'Haplotype_updated', 'rsID':'rsID_mapped'}, axis = 1)
data['rsID_mapped_sorted'] = data['rsID_mapped'].apply(lambda x: ','.join(sorted(x.split(','))))
data = data.drop(columns=['rsID_mapped'], axis=1)
mapped_df = pd.merge(madhu, data, on = 'rsID_mapped_sorted', how = 'inner', sort = False)
mapped_df.drop_duplicates(subset='Haplotype_updated', inplace=True)
mapped_df = mapped_df.drop(columns=['Zygosity'], axis=1)
mapped_df = mapped_df[['rsID_mapped_sorted', 'Zygosity_mapped', 'Haplotype_updated']]
print(mapped_df)