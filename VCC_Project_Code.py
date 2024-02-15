########################################### shell script for Extraction of the Required columns from the Dbsnp data ###################################
#zgrep "VC=SNV" dbSNP_2023_08_22.vcf.gz | awk 'BEGIN {OFS="\t"; print "CHROM", "POS", "rsID", "REF", "ALT"} {print $1, $2, $3, $4, $5}' > output_file.tsv


########################################### Extraction of Covered region variants from the Dbsnp data ###############################

import os
from tqdm import tqdm  # Import tqdm for progress bar

def read_bed_file(bed_file):
    bed_positions = set()
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                try:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    continue  # Skip this line if start or end position is not an integer
                for pos in range(start, end + 1):
                    bed_positions.add((chrom, pos))
    return bed_positions

def filter_tsv_file(tsv_file, bed_positions):
    filtered_tsv_records = []
    row_counter = 0  # Initialize a row counter
    with open(tsv_file, 'r') as f:
        for line in tqdm(f):  # Use tqdm to create a progress bar
            row_counter += 1  # Increment the row counter
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                try:
                    chrom = fields[0]
                    pos = int(fields[1])
                except ValueError:
                    continue  # Skip this line if 'POS' is not an integer
                if (chrom, pos) in bed_positions:
                    filtered_tsv_records.append(line)
    return filtered_tsv_records, row_counter  # Return the filtered records and row counter

def write_filtered_tsv(filtered_tsv_records, output_file, column_headers):
    with open(output_file, 'w') as f:
        f.write('\t'.join(column_headers) + '\n')
        for record in filtered_tsv_records:
            f.write(record)

def main():
    bed_file = r'C:/Users/GenepoweRx_Madhu/Downloads/BED_files/Covered_regions.bed'
    input_file = r'C:/Users/GenepoweRx_Madhu/Downloads/sen_spe_files_07_09_2023/output_file.tsv'
    output_file = r'C:/Users/GenepoweRx_Madhu/Downloads/sen_spe_files_07_09_2023/VCC_input_files/Truthset_dbsnp_SNV_covered_new.tsv'

    bed_positions = read_bed_file(bed_file)

    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))

    with open(input_file, 'r') as f:
        column_headers = f.readline().strip().split('\t')
    filtered_tsv_records, row_counter = filter_tsv_file(input_file, bed_positions)
    write_filtered_tsv(filtered_tsv_records, output_file, column_headers)

    print(f"Number of rows processed: {row_counter}")

if __name__ == "__main__":
    main()

####################################################### Covered regions variants extraction from the vcf files #######################

import os

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
    bed_file = r'C:/Users/GenepoweRx_Madhu/Downloads/BED_files/Covered_regions.bed'
    input_folder = r'C:/Users/GenepoweRx_Madhu/Downloads/sen_spe_files_07_09_2023/VCC_input_files/new_input_sample/new_bwa_files/'
    output_folder = r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/VCC_NEW/'

    bed_positions = read_bed_file(bed_file)

    vcf_files = [f for f in os.listdir(input_folder) if f.endswith('.vcf')]

    for vcf_file in vcf_files:
        input_file_path = os.path.join(input_folder, vcf_file)
        output_file_path = os.path.join(output_folder, vcf_file)
        filtered_vcf_records = filter_vcf_file(input_file_path, bed_positions)
        write_filtered_vcf(filtered_vcf_records, output_file_path)

if __name__ == "__main__":
    main()


################################### Comparision of the Truthset with the vcf files to get the required Table ##########################

import numpy as np
import pandas as pd
import os

# Function to calculate Jaccard similarity
def calculate_jaccard_similarity(set_a, set_b):
    intersection = len(set_a.intersection(set_b))
    union = len(set_a.union(set_b))
    return intersection / union if union != 0 else 0

# Define a function to process VCF files for a given tool and sample
def process_vcf(tool, sample, truth_set):
    vcf_path = f'/run/media/administrator/One Touch/New_vcc_data/VCC_NEW/{sample}_{tool}.vcf'
    
    if not os.path.exists(vcf_path):
        return None, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  # Return None and zeros for counts and metrics if the file does not exist
    
    vcf = pd.read_csv(vcf_path, comment='#', sep='\t', header=None, low_memory=False, encoding='latin-1')
    vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    
    def map_genotype(value):
        if '1/1' in value:
            return 'HOM'
        elif '0/1' in value:
            return 'HET'
        else:
            return 'HET'
    
    vcf['Zygosity'] = vcf['SAMPLE'].apply(map_genotype)
    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
    vcf['gnomADe_SAS_AF'] = vcf['CSQ'].str.split('|').str[56]
    vcf['gnomADe_AF'] = vcf['CSQ'].str.split('|').str[48]
    vcf['gnomADe_AF'] = vcf['gnomADe_AF'].replace('', 0).astype('float')
    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].replace('', 0)
    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].astype('float')
    vcf = vcf[~vcf['ALT'].str.contains(',')]
    #vcf = vcf[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO', 'FORMAT', 'SAMPLE', 'DP', 'VAF', 'gnomADe_AF', 'gnomADe_SAS_AF']]
    
    if tool == 'BCFTOOL':
        vcf['DP'] = vcf['INFO'].str.extract(r'DP=(\d+)')[0].fillna('0').astype(int)
    elif tool in ['DEEPVARIANT', 'HAPLOTYPECALLER']:
        vcf['DP'] = vcf['SAMPLE'].str.split(':').str[2].fillna('0').astype(int)
    elif tool in ['VARSCAN2', 'MUTECT2']:
        vcf['DP'] = vcf['SAMPLE'].str.split(':').str[3].fillna('0').astype(int)
        
    # Extract AD and calculate VAF based on the tool
    if tool == 'BCFTOOL':
        vcf['AD'] = vcf['SAMPLE'].str.split(':').str[6]
        vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
        vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
        vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
        vcf = vcf[vcf['QUAL'] >= 20]
    elif tool == 'VARSCAN2':
        vcf['RD'] = vcf['SAMPLE'].str.split(':').str[4].fillna('0').astype(int)
        vcf['AD'] = vcf['SAMPLE'].str.split(':').str[5].fillna('0').astype(int)
        vcf['VAF'] = vcf['AD'] / (vcf['RD'] + vcf['AD'])
    elif tool == 'HAPLOTYPECALLER':
        vcf['AD'] = vcf['SAMPLE'].str.split(':').str[1]
        vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
        vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
        vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
    elif tool == 'MUTECT2':
        vcf['AD'] = vcf['SAMPLE'].str.split(':').str[1]
        vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
        vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
        vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
    elif tool == 'DEEPVARIANT':
        vcf['AD'] = vcf['SAMPLE'].str.split(':').str[3]
        vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
        vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
        vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
        vcf = vcf[vcf['QUAL'] >= 20]
        
    before_count = len(vcf)
    vcf_vaf = vcf[vcf['VAF'] >= 0.1]
    vcf_after = vcf_vaf[vcf_vaf['DP'] >= 15]
    vcf_after = vcf_after[vcf_after['gnomADe_AF'] <= 0.5]
    vcf_after = vcf_after[vcf_after['gnomADe_SAS_AF'] <= 0.5]
    after_count = len(vcf_after)
    
    # Compare with the truth set
    truth_set['CHROM'] = truth_set['CHROM'].astype(str)
    truth_set['POS'] = truth_set['POS'].astype(int)
    truth_set['REF'] = truth_set['REF'].astype(str)
    truth_set['ALT'] = truth_set['ALT'].astype(str)

    truth_1 = truth_set.copy()
    truth_1['ALT'] = truth_1['ALT'].str.split(',')
    truth_1 = truth_1.explode('ALT')
    
    intersection_ab = pd.merge(vcf_after, truth_1, on=['CHROM', 'POS', 'REF', 'ALT'])
    intersection_ab = intersection_ab.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'])
    intersection_ab = len(intersection_ab)
    n_a = len(truth_set)
    n_b = len(vcf_after)
    
    TP = intersection_ab  # True Positives (A ∩ B)
    FP = n_b - TP  # False Positives (n(A) - TP)
    FN = n_a - TP  # False Negatives (n(B) - TP)
    TN = (1045445353 - 23761815) + FP 
    
    specificity = (TP) / (TP + TN + FP)
    sensitivity = TP / (TP + FN)
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    accuracy = (((TP + TN)) / (TP + TN + FP + FN))
    f1_score = 2 * ((precision * sensitivity) / (precision + sensitivity))
    
    # Calculate Jaccard Similarity
    truth_set_variants = set(zip(truth_set['CHROM'].astype(str), truth_set['POS'].astype(int), truth_set['REF'].astype(str), truth_set['ALT'].astype(str)))
    vcf_variants = set(zip(vcf_after['CHROM'].astype(str), vcf_after['POS'].astype(int), vcf_after['REF'].astype(str), vcf_after['ALT'].astype(str)))
    jaccard_similarity = calculate_jaccard_similarity(truth_set_variants, vcf_variants)
    
    # Calculate Concordance Rate
    concordance_rate = (TP) / (TP + FP + FN)
    
    return vcf, before_count, vcf_after, after_count, TP, TN, FP, FN, specificity, sensitivity, precision, recall, accuracy, f1_score, concordance_rate, jaccard_similarity

# Load the truth set in TSV format
truth_set_path = r'/run/media/administrator/One Touch/New_vcc_data/output.tsv'
truth_set = pd.read_csv(truth_set_path, sep='\t', header=0, low_memory=False, encoding='latin-1')
truth_set.columns = ['CHROM', 'POS', 'REF', 'ALT']

# Define a list of tools and samples
tools = ['BCFTOOL', 'VARSCAN2', 'MUTECT2', 'HAPLOTYPECALLER', 'DEEPVARIANT']
samples = ['12652705', '12652707', '12652709', '12652710', '12652712']

# Create a dictionary to store the results
results = {}

# Process each combination of tool and sample
for tool in tools:
    for sample in samples:
        vcf, before_count, vcf_after, after_count, TP, TN, FP, FN, specificity, sensitivity, precision, recall, accuracy, f1_score, concordance_rate, jaccard_similarity = process_vcf(tool, sample, truth_set)
        if vcf is not None:
            results[f'{sample}_{tool}'] = {
                'Sample': sample,
                'Tool': tool,
                'Before_Count': before_count,
                'After_Count(DP>=10)': after_count,
                'Before_After_Count_Difference': before_count - after_count,
                'HET_Count_Before': len(vcf[vcf['Zygosity'] == 'HET']),
                'HET_Count_After': len(vcf_after[vcf_after['Zygosity'] == 'HET']),
                'HET_Difference': len(vcf[vcf['Zygosity'] == 'HET']) - len(vcf_after[vcf_after['Zygosity'] == 'HET']),
                'HOM_Count_Before': len(vcf[vcf['Zygosity'] == 'HOM']),
                'HOM_Count_After': len(vcf_after[vcf_after['Zygosity'] == 'HOM']),
                'HOM_Difference': len(vcf[vcf['Zygosity'] == 'HOM']) - len(vcf_after[vcf_after['Zygosity'] == 'HOM']),
                'TP': TP,
                'TN': TN,
                'FP': FP,
                'FN': FN,
                'Specificity': specificity,
                'Sensitivity': sensitivity,
                'Precision': precision,
                'Recall': recall,
                'Accuracy': accuracy,
                'F1_Score': f1_score,
                'Concordance_Rate': concordance_rate,
                'Jaccard_Similarity': jaccard_similarity
            }

# Create a DataFrame from the results dictionary
results_df = pd.DataFrame.from_dict(results, orient='index')

# Insert the Truthset_count column at the 2nd index position
results_df.insert(2, 'Truthset_count', len(truth_set))
results_df = results_df.sort_values(by='Sample')

# Print the results DataFrame
results_df.to_excel(r'/run/media/administrator/One Touch/New_vcc_data/new_snp_0.1VAF_05_11_2023_BWA_DP15.xlsx', index=False) ### change the path
print(results_df)


################## True Positive rows mapping with ACMG columns for getting classification counts ###################################

truth_set_path = r'/run/media/administrator/One Touch/New_vcc_data/output.tsv'
truth_set = pd.read_csv(truth_set_path, sep='\t', header=0, low_memory=False, encoding='latin-1')
truth_1 = truth_set.copy()
truth_1['ALT'] = truth_1['ALT'].str.split(',')
truth_1 = truth_1.explode('ALT')

# Read your main data file
df = pd.read_excel(r'/run/media/administrator/One Touch/New_vcc_data/All_5samples.xlsx', engine='openpyxl') ### change the path
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
vaf_filters = [0.1]
dp_filters = [15]

# Create an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['Variant_classification', 'Tool', 'Sample', 'VAF_Filter', 'DP_Filter', 'Counts'])

# Iterate through combinations and calculate counts
for sample in samples:
    for tool in tools:
        for vaf_filter in vaf_filters:
            for dp_filter in dp_filters:
                # Read the corresponding VCF file for the tool
                vcf_file_path = f'/run/media/administrator/One Touch/New_vcc_data/VCC_NEW/{sample}_{tool}.vcf'
                vcf = pd.read_csv(vcf_file_path, comment='#', sep='\t', header=None, low_memory=False, encoding='latin-1')
                vcf.columns = ['CHROM', 'POS', 'rsID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

                # Extract VAF based on the tool
                if tool == 'BCFTOOL':
                    # Customize VAF extraction for BCFTOOL
                    vcf['DP'] = vcf['INFO'].str.extract(r'DP=(\d+)')[0].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[6]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
                    vcf['gnomADe_SAS_AF'] = vcf['CSQ'].str.split('|').str[56]
                    vcf['gnomADe_AF'] = vcf['CSQ'].str.split('|').str[48]
                    vcf['gnomADe_AF'] = vcf['gnomADe_AF'].replace('', 0).astype('float')
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].replace('', 0)
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].astype('float')
                    # For BCFTOOL, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                    filtered_data = filtered_data[filtered_data['gnomADe_AF'] <= 0.5]
                    filtered_data = filtered_data[filtered_data['gnomADe_SAS_AF'] <= 0.5]
                elif tool == 'VARSCAN2':
                    # Customize VAF extraction for VARSCAN2
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[3].fillna('0').astype(int)
                    vcf['RD'] = vcf['SAMPLE'].str.split(':').str[4].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[5].fillna('0').astype(int)
                    vcf['VAF'] = vcf['AD'] / (vcf['RD'] + vcf['AD'])
                    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
                    vcf['gnomADe_SAS_AF'] = vcf['CSQ'].str.split('|').str[56]
                    vcf['gnomADe_AF'] = vcf['CSQ'].str.split('|').str[48]
                    vcf['gnomADe_AF'] = vcf['gnomADe_AF'].replace('', 0).astype('float')
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].replace('', 0)
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].astype('float')
                    # Apply VAF and DP filters
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                    filtered_data = filtered_data[filtered_data['gnomADe_AF'] <= 0.5]
                    filtered_data = filtered_data[filtered_data['gnomADe_SAS_AF'] <= 0.5]
                elif tool == 'MUTECT2':
                    # Customize VAF extraction for MUTECT2
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[3].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[1]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
                    vcf['gnomADe_SAS_AF'] = vcf['CSQ'].str.split('|').str[56]
                    vcf['gnomADe_AF'] = vcf['CSQ'].str.split('|').str[48]
                    vcf['gnomADe_AF'] = vcf['gnomADe_AF'].replace('', 0).astype('float')
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].replace('', 0)
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].astype('float')
                    # For MUTECT2, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                    filtered_data = filtered_data[filtered_data['gnomADe_AF'] <= 0.5]
                    filtered_data = filtered_data[filtered_data['gnomADe_SAS_AF'] <= 0.5]
                elif tool == 'HAPLOTYPECALLER':
                    # Customize VAF extraction for HAPLOTYPECALLER
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[2].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[1]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
                    vcf['gnomADe_SAS_AF'] = vcf['CSQ'].str.split('|').str[56]
                    vcf['gnomADe_AF'] = vcf['CSQ'].str.split('|').str[48]
                    vcf['gnomADe_AF'] = vcf['gnomADe_AF'].replace('', 0).astype('float')
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].replace('', 0)
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].astype('float')
                    # For HAPLOTYPECALLER, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                    filtered_data = filtered_data[filtered_data['gnomADe_AF'] <= 0.5]
                    filtered_data = filtered_data[filtered_data['gnomADe_SAS_AF'] <= 0.5]
                elif tool == 'DEEPVARIANT':
                    # Customize VAF extraction for DEEPVARIANT
                    vcf['DP'] = vcf['SAMPLE'].str.split(':').str[2].fillna('0').astype(int)
                    vcf['AD'] = vcf['SAMPLE'].str.split(':').str[3]
                    vcf['RD'] = vcf['AD'].str.split(',').str[0].astype(int)
                    vcf['A_D'] = vcf['AD'].str.split(',').str[1].astype(int)
                    vcf['VAF'] = vcf['A_D'] / (vcf['RD'] + vcf['A_D'])
                    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
                    vcf['gnomADe_SAS_AF'] = vcf['CSQ'].str.split('|').str[56]
                    vcf['gnomADe_AF'] = vcf['CSQ'].str.split('|').str[48]
                    vcf['gnomADe_AF'] = vcf['gnomADe_AF'].replace('', 0).astype('float')
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].replace('', 0)
                    vcf['gnomADe_SAS_AF'] = vcf['gnomADe_SAS_AF'].astype('float')
                    # For DEEPVARIANT, you can choose not to filter based on 'DP'
                    filtered_data = vcf[vcf['VAF'] >= vaf_filter]
                    filtered_data = filtered_data[filtered_data['DP'] >= dp_filter]
                    filtered_data = filtered_data[filtered_data['gnomADe_AF'] <= 0.5]
                    filtered_data = filtered_data[filtered_data['gnomADe_SAS_AF'] <= 0.5]

                    
                filtered_data = pd.merge(filtered_data, truth_1, on = ['CHROM', 'POS', 'REF', 'ALT'])
                filtered_data = filtered_data.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'])

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

result_df.to_excel(r'/run/media/administrator/One Touch/New_vcc_data/BWA_MEM_True_positives_ACMG_05_11_2023.xlsx', index = False) ### change the path
# Display the final result
print(result_df)

###################################### Sensitivity, Specificity, Accuracy, Precision, F1_Score representation #########################

import matplotlib.pyplot as plt
import pandas as pd

data = {
    "Tool": ["HAPLOTYPECALLER", "MUTECT2", "BCFTOOL", "DEEPVARIANT", "VARSCAN2"],
    "Specificity": [0.0000265063251134942, 0.0000263955294386954, 0.0000251969847324033, 0.000025112424017798, 0.0000251108572198479]
}

df = pd.DataFrame(data)

# Define a list of colors to use for the data points
colors = ['blue', 'green', 'red', 'purple', 'orange']

# Create the line graph with lines and dots of different colors
for i in range(len(df) - 1):
    plt.plot(df["Tool"][i:i+2], df["Specificity"][i:i+2], marker='o', linestyle='-', color=colors[i])

# Add labels and title
plt.xlabel("Tool")
plt.ylabel("Specificity")
plt.title("")

# Rotate the x-axis labels for better readability
plt.xticks(rotation=0)

# Display the graph
plt.tight_layout()
#plt.savefig(r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/visualizations_new/BWA_0.1_DP15_Specificity_all_samples.tiff', dpi=600, bbox_inches='tight') ### change the path
plt.show()


############### Bowtie2, BWA_MEM data Sensitivity, Specificity, Accuracy, Precision, F1_Score representation ##########################

import matplotlib.pyplot as plt

# Data from the table
tools = ["VARSCAN2", "BCFTOOL", "HAPLOTYPECALLER", "MUTECT2", "DEEPVARIANT"]
sample_names = ["Sample 1"]

Specificity_BOWTIE2 = [
    [0.995393155, 0.995815396, 0.995907728, 0.994163425, 0.995865685]
]
Specificity_BWA_MEM = [
    [0.995337165, 0.994941576, 0.995894883, 0.992999941, 0.995939223]
]

# Create a single subplot for one sample
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Plot a line for specificity
ax.plot(range(len(tools)), Specificity_BOWTIE2[0], marker='o', label='BOWTIE2', color='b', alpha=0.6)

# Plot a line for sensitivity
ax.plot(range(len(tools)), Specificity_BWA_MEM[0], marker='o', label='BWA_MEM', color='r', alpha=0.6)

# Add labels and a legend
ax.set_xticks(range(len(tools)))
ax.set_xticklabels(tools)
ax.set_xlabel('Tools')
ax.set_ylabel('Precision')
ax.set_title("Precision")
ax.legend()

# Show the plot
plt.grid(True, axis='y', linestyle='--', alpha=0.6)
#plt.savefig(r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/visualizations_new/BWA_MEM_Bowtie2_Precision_all_samples.tiff', dpi=600, bbox_inches='tight') ### change the path
plt.show()

################### Individual Bowtie2, BWA_MEM data Sensitivity, Specificity, Accuracy, Precision, F1_Score representation ############

import matplotlib.pyplot as plt
import pandas as pd

# Data from the table
tools = ["12652705_VARSCAN2", "12652705_BCFTOOL", "12652705_HAPLOTYPECALLER", "12652705_MUTECT2", "12652705_DEEPVARIANT", "12652707_VARSCAN2", "12652707_BCFTOOL", "12652707_HAPLOTYPECALLER", "12652707_MUTECT2", "12652707_DEEPVARIANT", "12652709_VARSCAN2", "12652709_BCFTOOL", "12652709_HAPLOTYPECALLER", "12652709_MUTECT2", "12652709_DEEPVARIANT", "12652710_VARSCAN2", "12652710_BCFTOOL", "12652710_MUTECT2", "12652710_HAPLOTYPECALLER", "12652710_DEEPVARIANT", "12652712_VARSCAN2", "12652712_BCFTOOL", "12652712_HAPLOTYPECALLER", "12652712_MUTECT2", "12652712_DEEPVARIANT"]
Bowtie2_Sensitivity = [0.994506543, 0.99536965, 0.995241653, 0.994050902, 0.995584239, 0.995038122, 0.995443265, 0.995641806, 0.992692939, 0.995280289, 0.995555729, 0.996088861, 0.996010743, 0.994710039, 0.995857003, 0.995765298, 0.996004655, 0.996087331, 0.995084845, 0.996197558, 0.996100085, 0.996170549, 0.996557108, 0.994278398, 0.996409336]
BWA_MEM_Sensitivity = [0.994431094, 0.994852437, 0.995025779, 0.991482547, 0.995363629, 0.995162851, 0.994590696, 0.995595383, 0.991557009, 0.995666479, 0.995454016, 0.99523662, 0.996130116, 0.993670191, 0.99595625, 0.995624548, 0.995772395, 0.996385194, 0.994506884, 0.996179269, 0.996013315, 0.994255735, 0.996337945, 0.993783074, 0.996530485]

# Combine the data into a single DataFrame
data = {
    "Tool": tools,
    "Bowtie2_Sensitivity": Bowtie2_Sensitivity,
    "BWA_MEM_Sensitivity": BWA_MEM_Sensitivity
}
df = pd.DataFrame(data)

# Create a single subplot for the graph
fig, ax = plt.subplots(figsize=(12, 8))

# Plot lines for Bowtie2 Sensitivity
ax.plot(range(len(df)), df["Bowtie2_Sensitivity"], marker='o', label='Bowtie2', color='b', alpha=0.6)

# Plot lines for BWA_MEM Sensitivity
ax.plot(range(len(df)), df["BWA_MEM_Sensitivity"], marker='o', label='BWA_MEM', color='r', alpha=0.6)

# Set x-axis labels
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df["Tool"], rotation=90)

# Add labels and a legend
ax.set_xlabel('Tools')
ax.set_ylabel('Precision')
ax.set_title("Precision")
ax.legend()

# Show the plot
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()
#plt.savefig(r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/visualizations_new/Individual_Precision.png', dpi=600, bbox_inches='tight') ### change the path
plt.show()


############################### Overall Sensitivity vs Specificity representations #######################################################

import matplotlib.pyplot as plt
import pandas as pd

data = {
    "Tool": ["VARSCAN2", "BCFTOOL", "MUTECT2", "HAPLOTYPECALLER", "DEEPVARIANT"],
    "Specificity": [2.51109E-05, 2.5197E-05, 2.63955E-05, 2.65063E-05, 2.51124E-05],
    "Sensitivity": [0.001079716, 0.001083419, 0.001134955, 0.001139719, 0.001079783]
}

df = pd.DataFrame(data)

# Define a list of colors to use for the data points
colors = ['blue', 'green', 'red', 'purple', 'orange']

# Increase the figure dimensions by specifying figsize
plt.figure(figsize=(8, 5))

# Create a scatter plot with specificity on the x-axis and sensitivity on the y-axis
for i in range(len(df)):
    plt.scatter(df["Specificity"][i], df["Sensitivity"][i], color=colors[i], marker='o', label=df["Tool"][i])

# Add labels and title
plt.xlabel("Specificity")
plt.ylabel("Sensitivity")
plt.title("")

# Add a legend outside the plot
plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5))

# Display the graph
plt.grid(True)
plt.tight_layout()
#plt.savefig(r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/BWA_0.1_Specificity_vs_Sensitivity.tiff', dpi=600, bbox_inches='tight') ### change the path
plt.show()

############################# TPR vs FPR representation #####################################################################

import matplotlib.pyplot as plt

# Data for a single sample
sample_name = "All samples TPR vs FPR"
Tool = ['VARSCAN2', 'BCFTOOL', 'HAPLOTYPECALLER', 'MUTECT2', 'DEEPVARIANT']
TPR = [0.001079716, 0.001083419, 0.001139719, 0.001134955, 0.001079783]
FPR = [0.999974889, 0.999974803, 0.999973494, 0.999973604, 0.999974888]

# Create a single figure
plt.figure(figsize=(8, 6))

# Plot TPR vs. FPR for each tool with lines connecting the dots
for i in range(len(Tool)):
    plt.plot(FPR[i], TPR[i], marker='o', label=Tool[i], alpha=0.6)

# Add lines connecting the dots for each tool
for i in range(len(Tool) - 1):
    plt.plot([FPR[i], FPR[i+1]], [TPR[i], TPR[i+1]], linestyle='-', color='gray', alpha=0.6)

# Add labels and a legend to the right of the plot
plt.xlabel('False Positive Rate (FPR)')
plt.ylabel('True Positive Rate (TPR)')
plt.title(f'')
legend = plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))

# Show the plot
plt.grid(True, linestyle='--', alpha=0.6)
#plt.savefig(r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/BWA_0.1_ROC_TPR_FPR.tiff', dpi=600, bbox_inches='tight') ### change the path
plt.show()
