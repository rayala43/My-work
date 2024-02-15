import pandas as pd

# Function to read BED file and store genomic regions
def read_bed_file(bed_file):
    regions = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split('\t')
            start, end = int(start), int(end)
            if chrom not in regions:
                regions[chrom] = set()
            for pos in range(start, end + 1):
                regions[chrom].add(pos)
    return regions

# Read the BED file and store genomic regions
bed_file = r'C:/Users/GenepoweRx_Madhu/Downloads/BED_files/srinivas_sir_covered.bed'
regions = read_bed_file(bed_file)

# List of VCF files for father, mother, and son
vcf_files = ['12652705_BCFTOOL.vcf', '12652705_VARSCAN2.vcf', '12652705_HAPLOTYPECALLER.vcf', '12652705_MUTECT2.vaf', '12652705_DEEPVARIANT.vcf']

for vcf_file in vcf_files:
    # Create an empty list to store the matching variants for each family member
    matching_variants = []
    header_lines = []  # Store header lines from the VCF file

    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                header_lines.append(line)  # Capture header lines
                continue
            columns = line.strip().split('\t')
            chrom, position = columns[0], int(columns[1])
            if chrom in regions and position in regions[chrom]:
                matching_variants.append(columns)

    # Write the matching variants to an output VCF file
    output_filename = f'matching_positions{vcf_file.split("_")[2]}.vcf'
    
    # Write the header lines first
    with open(output_filename, 'w') as output_vcf:
        for header_line in header_lines:
            output_vcf.write(header_line)
    
    # Write the matching variants
    with open(output_filename, 'a') as output_vcf:
        for variant in matching_variants:
            output_vcf.write('\t'.join(variant) + '\n')