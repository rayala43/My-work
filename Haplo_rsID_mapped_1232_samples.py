import pandas as pd
import re
import sys
df = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/KHAIGPRX2_final_DP_Cond1.csv', sep = '\t')
df = df.rename({'rsID_y' : 'rsID'}, axis=1)
df['HET'] = df['INFO'].str.extract(r'HET=(\d)')
df['HOM'] = df['INFO'].str.extract(r'HOM=(\d)')
# Create a new column 'Zygosity' based on conditions
df['Zygosity'] = ''
df.loc[df['HOM'] == '1', 'Zygosity'] = 'Homozygous'
df.loc[df['HET'] == '1', 'Zygosity'] = 'Heterozygous'
df = df.drop(columns=['HET', 'HOM'], axis = 1)
data = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/Multiple_Positions.xlsx')
data = data.rename({'Haplotype':'Haplotype_updated', 'rsID':'rsID_mapped'}, axis = 1)
# Create new columns for mapped rsID and Zygosity
df['rsID_mapped'] = df.groupby('Haplotype')['rsID'].transform(lambda x: ','.join(set(x)))
# Create a Zygosity_mapped column where each row contains the corresponding Zygosity values for rsIDs in rsID_mapped
df['Zygosity_mapped'] = df.apply(lambda row: ','.join(row['Zygosity'] for rsID in row['rsID_mapped'].split(',')), axis=1)
# Create a new column with corresponding rsID and Zygosity pairs
df.drop_duplicates(subset='Haplotype', inplace=True)
mapped = pd.merge(df, data, on = 'rsID_mapped', how = 'inner', sort = False)
mapped.to_excel()
mapped