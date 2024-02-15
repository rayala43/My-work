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
    vcf_path = f'C:/Users/GenepoweRx_Madhu/Downloads/sen_spe_files_07_09_2023/VCC_input_files/{sample}_{tool}.vcf'
    
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
        vcf = vcf[vcf['QUAL'] >= 10]
        
    before_count = len(vcf)
    vcf_vaf = vcf[vcf['VAF'] >= 0.1]
    vcf_after = vcf_vaf[vcf_vaf['DP'] >= 18]
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
    
    intersection_ab = len(pd.merge(vcf_after, truth_1, on=['CHROM', 'POS', 'REF', 'ALT']))
    intersection_ab = intersection_ab.drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'])
    n_a = len(truth_set)
    n_b = len(vcf_after)
    
    TP = intersection_ab  # True Positives (A ∩ B)
    FP = n_b - TP  # False Positives (n(A) - TP)
    FN = n_a - TP  # False Negatives (n(B) - TP)
    TN = (1045445353 - 23761815) + FP 
    
    specificity = (TP) / (TP + TN + FN)
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
truth_set_path = r'C:/Users/GenepoweRx_Madhu/Downloads/COVERED_VCF_FILES_BED/output.tsv'
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
results_df.to_excel(r'C:/Users/GenepoweRx_Madhu/Downloads/sen_spe_files_07_09_2023/snp_0.1VAF_new_output_20_10_2023.xlsx', index=False)
print(results_df)
