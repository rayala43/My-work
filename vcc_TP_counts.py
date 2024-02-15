import numpy as np
import pandas as pd
import os
import sys
import re
import polars as pl
import xlsxwriter
import warnings
warnings.filterwarnings("ignore")
pd.set_option('display.max_columns', None)

truth_set_path = r'C:/Users/GenepoweRx_Madhu/Downloads/sen_spe_files_07_09_2023/VCC_input_files/Truthset_CHROM_POS_REF_ALT.tsv'
truth_set = pd.read_csv(truth_set_path, sep='\t', header=0, low_memory=False, encoding='latin-1')

# Read your main data file
df = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/MODY_new_05_10_2023/All_5samples.xlsx')
df['CHROM'] = 'chr' + df['Chr:Pos'].str.split(':').str[0]
df['POS'] = df['Chr:Pos'].str.split(':').str[1].astype(int)
df['REF'] = df['Ref/Alt'].str.split('/').str[0]
df['ALT'] = df['Ref/Alt'].str.split('/').str[1]
columns_to_check = ['CHROM', 'POS', 'REF', 'ALT']
# Use drop_duplicates to remove duplicated rows based on the specified columns
df = df.drop_duplicates(subset=columns_to_check)

# Define the sample, tool, VAF filter, and DP filter combinations
samples = ['12652705', '12652707', '12652709', '12652710', '12652712']
tools = ['BCFTOOL', 'VARSCAN2', 'MUTECT2', 'HAPLOTYPECALLER', 'DEEPVARIANT']
vaf_filters = [0.1, 0.01, 0.001]
dp_filters = [10, 18, 20, 22, 30]

# Create an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['Variant_classification', 'Tool', 'Sample', 'VAF_Filter', 'DP_Filter', 'Counts'])

# Iterate through combinations and calculate counts
for sample in samples:
    for tool in tools:
        for vaf_filter in vaf_filters:
            for dp_filter in dp_filters:
                # Read the corresponding VCF file for the tool
                vcf_file_path = f'C:/Users/GenepoweRx_Madhu/Downloads/sen_spe_files_07_09_2023/VCC_input_files/{sample}_{tool}.vcf'
                vcf = pd.read_csv(vcf_file_path, comment='#', sep='\t', header=None, low_memory=False, encoding='latin-1')
                vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

                # Extract VAF based on the tool
                if tool == 'BCFTOOL':
                    # Customize VAF extraction for BCFTOOL
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[6]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    # For BCFTOOL, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                elif tool == 'VARSCAN2':
                    # Customize VAF extraction for VARSCAN2
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[3].fillna('0').astype(int)
                    vcf['RD'] = vcf['SAMPLE'].str.split(':').str[4].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[5].fillna('0').astype(int)
                    vcf['VAF'] = vcf['AD'] / (vcf['RD'] + vcf['AD'])
                    # Apply VAF and DP filters
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                elif tool == 'MUTECT2':
                    # Customize VAF extraction for MUTECT2
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[3].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[1]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    # For MUTECT2, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                elif tool == 'HAPLOTYPECALLER':
                    # Customize VAF extraction for HAPLOTYPECALLER
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[2].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[1]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    # For HAPLOTYPECALLER, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                elif tool == 'DEEPVARIANT':
                    # Customize VAF extraction for DEEPVARIANT
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[2].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[3]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    # For DEEPVARIANT, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]

                    
                filtered_data = pd.merge(truth_set, filtered_data, on = ['CHROM', 'POS', 'REF', 'ALT'], how = 'inner', sort = False)

                # Merge the VCF data with the main data (df)
                merged_data = pd.merge(filtered_data, df, on=['CHROM', 'POS', 'REF', 'ALT'], how='left', sort=False)
                merged_data = merged_data.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'])

                # To fill NaN values in the 'Classification' column with 'NA', you need to assign the result back to the column.
                merged_data['Classification'] = merged_data['Classification'].fillna('NA')

                # List of your specified priorities
                priorities = ['Pathogenic', 'Likely Pathogenic', 'VUS/Weak Pathogenic', 'VUS/Conflicting', 'VUS', 'VUS/Weak Benign', 'Likely Benign', 'Benign']

                # Function to select the highest ranked value from Original_Column
                def select_highest_ranked(row):
                    values = row['Classification'].split(',')
                    for rank in priorities:
                        if rank in values:
                            return rank
                    return None  # Return None if no rank found

                # Apply the function to create the final selected rank column
                merged_data['Classification'] = merged_data.apply(select_highest_ranked, axis=1)
                
                # Calculate counts for each variant classification
                counts = merged_data['Classification'].value_counts(dropna=False)

                # Append the results to the result_df
                for classification, count in counts.items():
                    result_df = result_df.append({
                        'Variant_classification': classification,
                        'Tool': tool,
                        'Sample': sample,
                        'VAF_Filter': vaf_filter,
                        'DP_Filter': dp_filter,
                        'Counts': count
                    }, ignore_index=True)

# Reset index for the final result DataFrame
result_df.reset_index(drop=True, inplace=True)

result_df.to_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/True_positives_data_5samples_new_27_10_2023.xlsx', index = False)
# Display the final result
print(result_df)