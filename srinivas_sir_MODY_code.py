 ################################################# Covered TSV files from the bed pannel #################################################

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, MultipleLocator, ScalarFormatter


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
    with open(tsv_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                try:
                    chrom = fields[0]
                    pos = int(fields[1])
                except ValueError:
                    continue  # Skip this line if 'POS' is not an integer
                if (chrom, pos) in bed_positions:
                    filtered_tsv_records.append(line)
    return filtered_tsv_records

def write_filtered_tsv(filtered_tsv_records, output_file, column_headers):
    with open(output_file, 'w') as f:
        f.write('\t'.join(column_headers) + '\n')
        for record in filtered_tsv_records:
            f.write(record)

def main():
    bed_file = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files/Bed_file_covered.bed'
    input_folder = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files/'
    output_folder = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files/COVERED_TSV/Covered_new/'

    bed_positions = read_bed_file(bed_file)

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith('.tsv'):
            input_file = os.path.join(input_folder, filename)
            output_file = os.path.join(output_folder, filename)
            with open(input_file, 'r') as f:
                column_headers = f.readline().strip().split('\t')
            filtered_tsv_records = filter_tsv_file(input_file, bed_positions)
            write_filtered_tsv(filtered_tsv_records, output_file, column_headers)

if __name__ == "__main__":
    main()

################################# specific Matched Genes data for each sample from the 31 MODY samples ##################################

# Function to perform the matching and save the result as a TSV file
def process_file(input_path, output_path, df_gene, sample_name):
    vcf = pd.read_csv(input_path, sep="\t", low_memory=False)
    
    replacements = {'0/1': 'HET', '1/1': 'HOM'}
    vcf['Zygosity'] = vcf['GT'].replace(replacements)
    vcf = vcf[['ID', 'INFO', 'Zygosity']]
    
    vcf["Gene_Name"] = vcf["INFO"].str.extract('GENEINFO=(?P<GENEINFO>.+?);')
    vcf['Gene'] = vcf['Gene_Name'].apply(lambda x: ','.join([segment.split(':')[0] for segment in x.split('|')]) if pd.notnull(x) else '')
    
    vcf['Gene_Match'] = 'No'
    vcf['Matched_Gene'] = ''
    
    # Iterate through each gene in vcf['Gene']
    for index, genes in vcf['Gene'].iteritems():
        if isinstance(genes, str):
            gene_list = genes.split(',')
            for gene in gene_list:
                if gene in df_gene['Gene'].values:
                    vcf.at[index, 'Gene_Match'] = 'Yes'
                    vcf.at[index, 'Matched_Gene'] = gene
                    break
    
    matched = vcf[vcf['Gene_Match'] == 'Yes']
    
    # Add a new column with the sample name (extracted from the input file path)
    matched['Sample_Name'] = sample_name.split('_')[0]
    matched['CSQ'] = matched['INFO'].str.extract(r'CSQ=(.*)')
    matched = matched[['ID', 'INFO', 'CSQ', 'Zygosity', 'Matched_Gene', 'Sample_Name']]
    
    # Save the matched output as a TSV file in the specified output_path
    matched_output_filename = os.path.join(output_path, os.path.basename(input_path).replace('.tsv', '_matched_new.tsv'))
    matched.to_csv(matched_output_filename, sep="\t", index=False)

# Load the genes from the Excel file
df_gene = pd.read_excel(r'C:/Users/GenepoweRx_Madhu/Desktop/MODY_genes.xlsx')

# Set the input and output folder paths
input_folder = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files_31_mody/COVERED_TSV/Covered_new/'
output_folder = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files_31_mody/COVERED_TSV/Covered_new/matched/'

# Iterate through each TSV file in the input folder and process it
for filename in os.listdir(input_folder):
    if filename.endswith('.tsv'):
        input_file_path = os.path.join(input_folder, filename)
        sample_name = os.path.splitext(filename)[0]  # Extract the file name without extension
        process_file(input_file_path, output_folder, df_gene, sample_name)

print("Matching and saving process completed.")

############################################### Concatenate the 31 MODY samples Covered Genes ###########################################

def concatenate_tsv_files(folder_path, output_file):
    # Get a list of all TSV files in the folder
    tsv_files = [file for file in os.listdir(folder_path) if file.endswith('.tsv')]

    # Read each TSV file into a DataFrame and concatenate them
    dataframes = [pd.read_csv(os.path.join(folder_path, file), delimiter='\t') for file in tsv_files]
    concatenated_df = pd.concat(dataframes)

    # Save the concatenated DataFrame to a new TSV file
    concatenated_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    folder_path = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files_31_mody/COVERED_TSV/Covered_new/matched/'  # Replace this with the actual path to your TSV files
    output_file = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files_31_mody/COVERED_TSV/Covered_new/matched/output_concatenated.tsv'  # Replace this with the desired output file name

    concatenate_tsv_files(folder_path, output_file)


### CLIN_SIG, Consequence, IMPACT, EXON, INTRON, HGVSC, HGVSP, cDNA_position, Protein_position, Amino_acids, Codons columns extraction ###
####### Generate the selective Gene wise variants and columns extraction and save the files for each gene with HOM and HET sheets ########

df = pd.read_csv(r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files_31_mody/COVERED_TSV/Covered_new/matched/output_concatenated.tsv', sep='\t')
genes_to_drop = ['HLA-DQB1', 'HLA-DRB1', 'HLA-DQA1']
data = df[~df['Matched_Gene'].isin(genes_to_drop)]
df_new = data[['Sample_Name', 'ID', 'Zygosity', 'Matched_Gene', 'CSQ']]
df_new['gnomADg_AF'] = df_new['CSQ'].str.split('|').str[57]
df_new['gnomADg_SAS_AF'] = df_new['CSQ'].str.split('|').str[67]
df_new['csq'] = df_new['CSQ'].str.split(',')
df_new = df_new.explode('csq')
df_new['CLIN_SIG'] = df_new['csq'].str.split('|').str[70]
df_new['Consequence'] = df_new['csq'].str.split('|').str[1]
df_new['IMPACT'] = df_new['csq'].str.split('|').str[2]
df_new['EXON'] = df_new['csq'].str.split('|').str[8]
df_new['INTRON'] = df_new['csq'].str.split('|').str[9]
df_new['HGVSC'] = df_new['csq'].str.split('|').str[10]
df_new['HGVSP'] = df_new['csq'].str.split('|').str[11]
df_new['cDNA_position'] = df_new['csq'].str.split('|').str[12]
df_new['Protein_position'] = df_new['csq'].str.split('|').str[14]
df_new['Amino_acids'] = df_new['csq'].str.split('|').str[15]
df_new['Codons'] = df_new['csq'].str.split('|').str[16]
impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
# Sort the DataFrame based on the priority
df_new['IMPACT'] = pd.Categorical(df_new['IMPACT'], categories=impact_order, ordered=True)
df_new.sort_values(['Sample_Name', 'ID', 'Matched_Gene', 'IMPACT'], ascending=[True, True, True, True], inplace=True)

# Drop duplicates while keeping the first occurrence
filtered_df = df_new.drop_duplicates(subset=['Sample_Name', 'ID', 'Matched_Gene'], keep='first')
Cutoff_data = filtered_df[filtered_df['gnomADg_AF'] >= 0.02]
# Get unique gene classes
unique_genes = yes_new['Matched_Gene'].unique()
# Define the output directory
output_directory = r'C:\Users\GenepoweRx_Madhu\Downloads'

# Iterate through each unique gene and create separate Excel files
for gene in unique_genes:
    gene_data = yes_new[yes_new['Matched_Gene'] == gene]
    het_data = gene_data[gene_data['Zygosity'] == 'HET']
    hom_data = gene_data[gene_data['Zygosity'] == 'HOM']
    
    # Create a Pandas Excel writer
    output_filename = os.path.join(output_directory, f'{gene}_output_data.xlsx')
    writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')
    
    # Write each DataFrame to a separate Excel sheet
    het_data.to_excel(writer, sheet_name='HET', index=False)
    hom_data.to_excel(writer, sheet_name='HOM', index=False)
    
    # Save the Excel file
    writer.save()
    
    print(f"Excel file for {gene} created successfully at {output_filename}")


############################## All genes HET and HOM count (gnomADg_AF == Global, 0.001, 0.02) ###########################################


# Folder path containing the TSV files
folder_path = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files/COVERED_TSV/Covered_new/'

# Get a list of all files in the folder
tsv_files = [file for file in os.listdir(folder_path) if file.endswith('.tsv')]

# Create a list to store the data for the Excel file
excel_data = []

# Loop through each TSV file and process them
for tsv_file in tsv_files:
    # Read the TSV file into a DataFrame
    vcf = pd.read_csv(os.path.join(folder_path, tsv_file), sep="\t", low_memory=False)
    
    # Extract the 'CSQ' column
    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
    vcf['HET'] = vcf['INFO'].str.extract(r'HET=(\d)')
    vcf['HOM'] = vcf['INFO'].str.extract(r'HOM=(\d)')
    
    # Split the 'CSQ' column to get the desired columns
    csq_values = vcf['CSQ'].str.split('|')
    vcf['gnomADg_AF'] = csq_values.str[57].replace('', 0).astype(float)
    vcf['gnomADg_SAS_AF'] = csq_values.str[67].replace('', 0).astype(float)
    
    # Select only the desired columns
    vcf = vcf[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT',
               'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RDF', 'RDR', 'ADF',
               'ADR', 'gnomADg_AF', 'gnomADg_SAS_AF', 'HET', 'HOM']]
    
    # Filter rows with gnomADg_AF >= 0.02 (or) 0.001 as per our required range.
    vcf = vcf[vcf['gnomADg_AF'] >= 0.02]
    
    # Get the filename without the extension
    filename = os.path.splitext(tsv_file)[0]
    
    # Get the number of rows in the data
    num_rows = vcf.shape[0]
    
    # Count occurrences of 'HET = 1' and 'HOM = 1'
    het_count = vcf[vcf['HET'] == '1'].shape[0]
    hom_count = vcf[vcf['HOM'] == '1'].shape[0]
    
    # Append the data to the Excel data list
    excel_data.append([filename, num_rows, het_count, hom_count])

# Create a DataFrame from the Excel data list
excel_df = pd.DataFrame(excel_data, columns=['Filename', 'NumRows', 'HET=1_Count', 'HOM=1_Count'])

# Save the DataFrame as an Excel file
excel_file_path = os.path.join(folder_path, 'summary_0.02_Covered_positions.xlsx')
excel_df.to_excel(excel_file_path, index=False)

############################# Selective genes HET and HOM count (gnomADg_AF == Global, 0.001, 0.02) ######################################

# Folder path containing the TSV files
folder_path = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files/'

# Get a list of all files in the folder
tsv_files = [file for file in os.listdir(folder_path) if file.endswith('.tsv')]

# Create a list to store the data for the Excel file
excel_data = []

# Genes to consider for filtering
#target_genes = ['HNF1A', 'CEL', 'HNF1B', 'PAX4', 'HNF4A', 'INS', 'GCK', 'BLK', 'NEUROD1',
#                'ABCC8', 'PDX1', 'KCNJ11', 'KLF11', 'APPL1']
target_genes = ['HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1']

# Loop through each TSV file and process them
for tsv_file in tsv_files:
    # Read the TSV file into a DataFrame
    vcf = pd.read_csv(os.path.join(folder_path, tsv_file), sep="\t", low_memory=False)
    
    # Extract the 'CSQ' column
    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
    
    # Split the 'CSQ' column to get the desired columns
    csq_values = vcf['CSQ'].str.split('|')
    vcf['gnomADg_AF'] = csq_values.str[57].replace('', 0).astype(float)
    vcf['gnomADg_SAS_AF'] = csq_values.str[67].replace('', 0).astype(float)
    
    # Select only the desired columns
    vcf = vcf[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT',
               'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RDF', 'RDR', 'ADF',
               'ADR', 'gnomADg_AF', 'gnomADg_SAS_AF']]
    
    # Filter rows with gnomADg_AF >= 0.001
    #vcf = vcf[vcf['gnomADg_AF'] >= 0.001]
    
    # Extract gene names and create a new column 'Gene Name'
    vcf["Gene_Name"] = vcf["INFO"].str.extract('GENEINFO=(?P<GENEINFO>.+?);')
    vcf['Gene Name'] = vcf['Gene_Name'].apply(lambda x: ','.join([segment.split(':')[0] for segment in x.split('|')]) if pd.notnull(x) else '')
    
    # Filter rows containing target genes
    target_gene_rows = vcf[vcf['Gene Name'].isin(target_genes)]
    
    # Calculate gene counts and store in a dictionary
    gene_counts = target_gene_rows['Gene Name'].value_counts().to_dict()
    
    # Calculate HET and HOM counts for each gene
    het_hom_counts = {}
    for gene in target_genes:
        het_count = target_gene_rows.loc[target_gene_rows['Gene Name'] == gene, 'GT'].str.contains('1/0|0/1').sum()
        hom_count = target_gene_rows.loc[target_gene_rows['Gene Name'] == gene, 'GT'].str.contains('1/1').sum()
        het_hom_counts[gene] = (het_count, hom_count)
    
    # Get the filename without the extension
    filename = os.path.splitext(tsv_file)[0]
    
    # Append the data to the Excel data list
    for gene, count in gene_counts.items():
        het_count, hom_count = het_hom_counts[gene]
        excel_data.append([filename, gene, count, het_count, hom_count])

# Create a DataFrame from the Excel data list
excel_df = pd.DataFrame(excel_data, columns=['Filename', 'Gene Name', 'Gene Count', 'HET=1_Count', 'HOM=1_Count'])

# Save the DataFrame as an Excel file
excel_file_path = os.path.join(folder_path, 'summary_3genes_Global.xlsx')
excel_df.to_excel(excel_file_path, index=False)

####################################### Selective genes Clin_sig Extraction from the CSQ of the INFO #####################################################

# Folder path containing the TSV files
folder_path = r'C:/Users/GenepoweRx_Madhu/Downloads/Depth_files/'

# Get a list of all files in the folder
tsv_files = [file for file in os.listdir(folder_path) if file.endswith('.tsv')]

# Create a list to store the data for the Excel file
excel_data = []

# Genes to consider for filtering
target_genes = ['HNF1A', 'CEL', 'HNF1B', 'PAX4', 'HNF4A', 'INS', 'GCK', 'BLK', 'NEUROD1',
                'ABCC8', 'PDX1', 'KCNJ11', 'KLF11', 'APPL1']

# Loop through each TSV file and process them
for tsv_file in tsv_files:
    # Read the TSV file into a DataFrame
    vcf = pd.read_csv(os.path.join(folder_path, tsv_file), sep="\t", low_memory=False)
    
    # Extract the 'CSQ' column
    vcf['CSQ'] = vcf['INFO'].str.extract(r'CSQ=(.*)')
    
    # Split the 'CSQ' column to get the desired columns
    csq_values = vcf['CSQ'].str.split(',')
    vcf['csq'] = csq_values.str[0]
    vcf['CLIN_SIG'] = vcf['csq'].str.split('|').str[70]
    vcf['CLIN_SIG_new'] = vcf['CLIN_SIG'].str.split('&').str[0]
    vcf['CLIN_SIG_new_Final'] = vcf['CLIN_SIG_new'].str.split('/').str[0]
    
    # Split the 'CSQ' column to get the desired columns
    csq_values = vcf['CSQ'].str.split('|')
    vcf['gnomADg_AF'] = csq_values.str[57].replace('', 0).astype(float)
    vcf['gnomADg_SAS_AF'] = csq_values.str[67].replace('', 0).astype(float)
    
    # Select only the desired columns
    vcf = vcf[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'GT',
               'GQ', 'SDP', 'DP', 'RD', 'AD', 'FREQ', 'PVAL', 'RDF', 'RDR', 'ADF',
               'ADR', 'gnomADg_AF', 'gnomADg_SAS_AF', 'CSQ', 'CLIN_SIG_new_Final']]
    
    # Filter rows with gnomADg_AF >= 0.02
    vcf = vcf[vcf['gnomADg_AF'] >= 0.02]
    
    # Extract gene names and create a new column 'Gene Name'
    vcf["Gene_Name"] = vcf["INFO"].str.extract('GENEINFO=(?P<GENEINFO>.+?);')
    vcf['Gene Name'] = vcf['Gene_Name'].apply(lambda x: ','.join([segment.split(':')[0] for segment in x.split('|')]) if pd.notnull(x) else '')
    
    # Filter rows containing target genes
    target_gene_rows = vcf[vcf['Gene Name'].isin(target_genes)]
    
    # Calculate gene counts and store in a dictionary
    gene_counts = target_gene_rows['Gene Name'].value_counts().to_dict()
    
    # Calculate HET and HOM counts for each gene
    het_hom_counts = {}
    het_clin_sig_counts = {}
    hom_clin_sig_counts = {}
    for gene in target_genes:
        het_count = target_gene_rows.loc[target_gene_rows['Gene Name'] == gene, 'GT'].str.contains('1/0|0/1').sum()
        hom_count = target_gene_rows.loc[target_gene_rows['Gene Name'] == gene, 'GT'].str.contains('1/1').sum()
        het_hom_counts[gene] = (het_count, hom_count)
        het_clin_sig_counts[gene] = target_gene_rows.loc[(target_gene_rows['Gene Name'] == gene) & (target_gene_rows['GT'].str.contains('1/0|0/1')), 'CLIN_SIG_new_Final'].value_counts().to_dict()
        hom_clin_sig_counts[gene] = target_gene_rows.loc[(target_gene_rows['Gene Name'] == gene) & (target_gene_rows['GT'].str.contains('1/1')), 'CLIN_SIG_new_Final'].value_counts().to_dict()
    
    # Get the filename without the extension
    filename = os.path.splitext(tsv_file)[0]
    
    # Append the data to the Excel data list
    for gene, count in gene_counts.items():
        het_count, hom_count = het_hom_counts[gene]
        het_clin_sig_counts_gene = het_clin_sig_counts.get(gene, {})
        hom_clin_sig_counts_gene = hom_clin_sig_counts.get(gene, {})
        for clin_sig in het_clin_sig_counts_gene.keys():
            het_clin_count = het_clin_sig_counts_gene.get(clin_sig, 0)
            hom_clin_count = hom_clin_sig_counts_gene.get(clin_sig, 0)
            excel_data.append([filename, gene, count, het_count, hom_count, clin_sig, het_clin_count, hom_clin_count])

# Create a DataFrame from the Excel data list
excel_df = pd.DataFrame(excel_data, columns=['Filename', 'Gene Name', 'Gene Count', 'HET=1_Count', 'HOM=1_Count', 'CLIN_SIG', 'HET_CLIN_Count', 'HOM_CLIN_Count'])

# Save the DataFrame as an Excel file
excel_file_path = os.path.join(folder_path, 'summary_new_clin_sig.xlsx')
excel_df.to_excel(excel_file_path, index=False)


################################################# HNF1A Gene Boxplot for gnomADg_AF and gnomADg_SAS_AF by Zygosity ######################

# Assuming you have loaded your data into the 'filtered_data' DataFrame

# Group the data by 'Zygosity'
grouped_data = filtered_data.groupby('Zygosity')

# Create a larger figure with higher DPI
fig = plt.figure(figsize=(3, 4), dpi=2000)

# Create positions for HET and HOM boxplots
het_positions = [0.5, 0.8]
hom_positions = [1.5, 1.8]

# Extract the data for HET and HOM
het_gnomADg_AF = grouped_data.get_group('HET')['gnomADg_AF']
het_gnomADg_SAS_AF = grouped_data.get_group('HET')['gnomADg_SAS_AF']
hom_gnomADg_AF = grouped_data.get_group('HOM')['gnomADg_AF']
hom_gnomADg_SAS_AF = grouped_data.get_group('HOM')['gnomADg_SAS_AF']

# Create boxplots for HET and HOM with reduced width
boxplot = plt.boxplot([het_gnomADg_AF, het_gnomADg_SAS_AF, hom_gnomADg_AF, hom_gnomADg_SAS_AF],
                      positions=[0.5, 0.8, 1.5, 1.8],
                      labels=['HET_AF', 'HET_SAS_AF', 'HOM_AF', 'HOM_SAS_AF'],
                      patch_artist=True,  # Enable patch_artist for colored boxes
                      widths=0.15)  # Adjust the width of the boxplots

# Set box colors
colors = ['#FF6103', '#FF6103', '#1E90FF', '#1E90FF']  # Distinct colors for HET and HOM
for box, color in zip(boxplot['boxes'], colors):
    box.set(facecolor=color)

# Set x-axis labels and title
plt.title('HNF1A Gene Boxplot for gnomADg_AF and gnomADg_SAS_AF by Zygosity', fontsize=4.8)
plt.ylabel('Value', fontsize=5)
plt.xlabel('Zygosity', fontsize=5)

# Set header text size for x-axis labels
plt.xticks(fontsize=3.8)
plt.yticks(fontsize=5.5)

# Adjust x-axis limits to decrease space
plt.xlim(0.2, 2.1)

# Calculate medians
medians = [median.get_ydata().mean() for median in boxplot['medians']]

# Adding median values as text
#for pos, median in zip([1, 2, 4, 5], medians):
#    plt.text(pos, median, f"{median:.2f}", color='black', ha='center', va='bottom', fontsize=8)

# Save the figure with higher quality
#plt.savefig(r'C:\Users\GenepoweRx_Madhu\Downloads\Depth_files_31_mody\COVERED_TSV\Covered_new\matched\14_genes_new_files\HET_HOM_boxplot_new.png', dpi=4000)

# Show the plot
plt.show()

############################# 12652713 and 17751406 samples Zygosity wise percentile breaken axis Bar graph ##############################

#df = Percentile data from the analysis

# Create positions for bars
x = np.arange(len(df))

# Set width of bars
width = 0.35

# Create the figure and axis
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [100, 1]})
plt.subplots_adjust(hspace=0.05)

# Bar chart for Hetero_count (make bars for certain categories wider)
hetero_bars = ax1.bar(x - width/2, df['Hetero_count'], width, label='Hetero Count')

# Bar chart for Homo_count (overlaying on ax1)
homo_bars = ax1.bar(x + width/2, df['Homo_count'], width, label='Homo Count')

# Apply different y-axis scaling functions for each range
def scale_function(y, category):
    if category == "(>25-50)":
        return y / 0.09
    elif category in ["(>50-75)", ">75"]:
        return y / 0.018
    else:
        return y

ax1.yaxis.set_major_locator(LogLocator(subs=[1.0]))
ax1.yaxis.set_major_formatter(ScalarFormatter())
ax1.yaxis.set_minor_locator(MultipleLocator(base=1))

for i, bar in enumerate(hetero_bars):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2, height, str(height), ha='center', va='bottom', fontsize=8)

for i, bar in enumerate(homo_bars):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2, height, str(height), ha='center', va='bottom', fontsize=8)

# Apply scaling functions to the counts based on categories
for i, (h_bar, g_bar) in enumerate(zip(hetero_bars, homo_bars)):
    h_height = h_bar.get_height()
    g_height = g_bar.get_height()
    category = df['Percentile'][i]
    h_bar.set_height(scale_function(h_height, category))
    g_bar.set_height(scale_function(g_height, category))

# Setting y-axis label for the top subplot
ax1.set_ylabel('Count')

# Create a broken y-axis for the bottom subplot
ax2.spines['top'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.tick_params(labeltop=False)

d = .015
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot([-d, +d], [-d, +d], **kwargs)
ax1.plot([1 - d, 1 + d], [-d, +d], **kwargs)

kwargs.update(transform=ax2.transAxes)
ax2.plot([-d, +d], [1 - d, 1 + d], **kwargs)
ax2.plot([1 - d, 1 + d], [1 - d, 1 + d], **kwargs)

# Setting y-axis label for the bottom subplot
ax2.set_ylabel('')

# Set x-axis ticks and labels
ax2.set_xticks(x)
ax2.set_xticklabels(df['Percentile'])

# Set custom x-axis labels
custom_labels = ["(0-25)", "(>25-50)", "(>50-75)", ">75"]
ax2.set_xticklabels(custom_labels)

ax2.set_xlabel('Percentile')

# Adding legend to the top subplot
ax1.legend()

plt.savefig(r'C:/Users/GenepoweRx_Madhu/Downloads/Final_17751406.svg', dpi=250, bbox_inches='tight')

# Show the plot
plt.show()

############################################### 12652713 and 17751406 samples Zygosity wise percentile bargraph #########################

# Create positions for bars
x = np.arange(len(df))

# Set width of bars
width = 0.35

# Create the figure and axis
fig, ax = plt.subplots()

# Create bar chart for Hetero_count
bar1 = ax.bar(x - width/2, df['Hetero_count'], width, label='Hetero Count')

# Create bar chart for Homo_count
bar2 = ax.bar(x + width/2, df['Homo_count'], width, label='Homo Count')

# Set x-axis ticks and labels
ax.set_xticks(x)
ax.set_xticklabels(df['Percentile'])
ax.set_xlabel('Percentile')

# Set y-axis label
ax.set_ylabel('Count')

# Adding data labels
def add_labels(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

add_labels(bar1)
add_labels(bar2)

# Adding legend
ax.legend()

#plt.savefig(r'C:/Users/GenepoweRx_Madhu/Downloads/3_genes_new_17751406.svg', dpi=450, bbox_inches='tight')

# Show the plot
plt.show()

############################### Extracting the Functional aspects columns from VEP and Group them by variant, Location and ALT ###########

df = pd.read_excel(r'C:/Variant data collect from the VEP for the specific genes and InDel Files.xlsx')
data = pd.read_excel(r'C:/variant and Location and Allele data from the all 31 mody samples.xlsx')
merged_df = pd.merge(df, data, on = ['Uploaded_variation', 'Location'], how = 'left', sort=False)
matching_rows = merged_df[merged_df['Allele'] == merged_df['ALT']]
grouped_data = matching_rows.groupby(['Uploaded_variation', 'Location', 'Allele']).agg(lambda x: ', '.join(x.unique())).reset_index()
