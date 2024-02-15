import pandas as pd

# Define the genes to check
genes_to_check = ["ABCG2", "CACNA1S", "CFTR", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "DPYD", "G6PD", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "OR4F5"]

# Initialize an empty list to store data
gene_info_list = []

# Open the VCF file
with open('C:/Users/GenepoweRx_Madhu/Downloads/718_final.vcf', 'r') as vcf_file:
    lines = vcf_file.readlines()

    # Initialize a record flag
    is_in_record = False

    # Initialize a dictionary to store record data
    record_data = {}

    for line in lines:
        line = line.strip()

        if line.startswith('#'):
            continue  # Skip header lines

        if line == "":
            if is_in_record:
                gene_symbols = record_data.get("INFO", {}).get("GENEINFO", "").split('|')
                for gene_symbol in gene_symbols:
                    gene_data = gene_symbol.split(':')
                    gene_name = gene_data[0]

                    if gene_name in genes_to_check:
                        # Extract relevant fields and store them as a dictionary
                        data = {
                            "Gene": gene_name,
                            "#CHROM": record_data.get("#CHROM", ""),
                            "POS": record_data.get("POS", ""),
                            "ID": record_data.get("ID", ""),
                            "REF": record_data.get("REF", ""),
                            "ALT": record_data.get("ALT", ""),
                            "QUAL": record_data.get("QUAL", ""),
                            "FILTER": record_data.get("FILTER", ""),
                            "INFO": record_data.get("INFO", ""),
                            "FORMAT": record_data.get("FORMAT", ""),
                            "Sample1": record_data.get("Sample1", ""),
                            "GT": record_data.get("GT", ""),
                            "GQ": record_data.get("GQ", ""),
                            "SDP": record_data.get("SDP", ""),
                            "DP": record_data.get("DP", ""),
                            "RD": record_data.get("RD", ""),
                            "AD": record_data.get("AD", ""),
                            "FREQ": record_data.get("FREQ", ""),
                            "PVAL": record_data.get("PVAL", ""),
                            "RBQ": record_data.get("RBQ", ""),
                            "ABQ": record_data.get("ABQ", ""),
                            "RDF": record_data.get("RDF", ""),
                            "RDR": record_data.get("RDR", ""),
                            "ADF": record_data.get("ADF", ""),
                            "ADR": record_data.get("ADR", ""),
                        }

                        # Extract CSQ if it's present in INFO
                        if "CSQ" in record_data.get("INFO", {}):
                            csq_values_list = record_data["INFO"]["CSQ"].split(',')
                            csq_values = csq_values_list[0].split('|')  # Corrected this line

                            csq_columns = [
                                "Allele",
                                "Consequence",
                                "Impact",
                                "Symbol",
                                "Gene",
                                "Feature_type",
                                "Feature",
                                "Biotype",  # Adding Biotype here
                                "EXON",
                                "INTRON",
                                "HGVSc",
                                "HGVSp",
                                "cDNA_position",
                                "CDS_position",
                                "Protein_position",
                                "Amino_acids",
                                "Codons",
                                "Existing_variation",
                                "DISTANCE",
                                "STRAND",
                                "FLAGS",
                                "SYMBOL_SOURCE",
                                "HGNC_ID",
                                "CANONICAL",
                                "TSL",
                                "APPRIS",
                                "CCDS",
                                "ENSP",
                                "SWISSPROT",
                                "TREMBL",
                                "UNIPARC",
                                "REFSEQ_MATCH",
                                "GENE_PHENO",
                                "SIFT",
                                "PolyPhen",
                                "DOMAINS",
                                "HGVS_OFFSET",
                                "GMAF",
                                "AFR_MAF",
                                "AMR_MAF",
                                "EAS_MAF",
                                "EUR_MAF",
                                "SAS_MAF",
                                "AA_MAF",
                                "EA_MAF",
                                "gnomAD_AF",
                                "gnomAD_AFR_AF",
                                "gnomAD_AMR_AF",
                                "gnomAD_ASJ_AF",
                                "gnomAD_EAS_AF",
                                "gnomAD_FIN_AF",
                                "gnomAD_NFE_AF",
                                "gnomAD_OTH_AF",
                                "gnomAD_SAS_AF",
                                "MAX_AF",
                                "MAX_AF_POPS",
                                "CLIN_SIG",
                                "SOMATIC",
                                "PHENO",
                                "PUBMED",
                                "MOTIF_NAME",
                                "MOTIF_POS",
                                "HIGH_INF_POS",
                                "MOTIF_SCORE_CHANGE",
                                "TRANSCRIPT_ID"
                            ]

                            csq_data = {column: value if value != '' else 'NA' for column, value in
                                        zip(csq_columns, csq_values)}
                            data.update({"CSQ_" + column: csq_data.get(column, "NA") for column in csq_columns})

                            # Extract the "Biotype" values and concatenate them
                            biotype_values = []
                            for transcript_csq in csq_values_list:
                                biotype = transcript_csq.split('|')[7]
                                biotype_values.append(biotype)

                            # Concatenate the Biotype values with commas
                            data["Biotype"] = ",".join(biotype_values)

                            gene_info_list.append(data)

                # Reset the record_data dictionary and flag
                record_data = {}
                is_in_record = False
        else:
            # Split the line into fields
            fields = line.split('\t')
            if len(fields) >= 10:
                # Store the fields in the record_data dictionary
                record_data["#CHROM"] = fields[0]
                record_data["POS"] = fields[1]
                record_data["ID"] = fields[2]
                record_data["REF"] = fields[3]
                record_data["ALT"] = fields[4]
                record_data["QUAL"] = fields[5]
                record_data["FILTER"] = fields[6]
                record_data["INFO"] = fields[7]
                record_data["FORMAT"] = fields[8]
                record_data["Sample1"] = fields[9]

                # Extract GT, GQ, SDP, DP, RD, AD, FREQ, PVAL, RBQ, ABQ, RDF, RDR, ADF, ADR values
                format_fields = fields[8].split(':')
                sample_fields = fields[9].split(':')
                record_data["GT"] = sample_fields[format_fields.index("GT")]
                record_data["GQ"] = sample_fields[format_fields.index("GQ")]
                record_data["SDP"] = sample_fields[format_fields.index("SDP")]
                record_data["DP"] = sample_fields[format_fields.index("DP")]
                record_data["RD"] = sample_fields[format_fields.index("RD")]
                record_data["AD"] = sample_fields[format_fields.index("AD")]
                record_data["FREQ"] = sample_fields[format_fields.index("FREQ")]
                record_data["PVAL"] = sample_fields[format_fields.index("PVAL")]
                record_data["RBQ"] = sample_fields[format_fields.index("RBQ")]
                record_data["ABQ"] = sample_fields[format_fields.index("ABQ")]
                record_data["RDF"] = sample_fields[format_fields.index("RDF")]
                record_data["RDR"] = sample_fields[format_fields.index("RDR")]
                record_data["ADF"] = sample_fields[format_fields.index("ADF")]
                record_data["ADR"] = sample_fields[format_fields.index("ADR")]

                is_in_record = True

# Create a DataFrame from the list
gene_info = pd.DataFrame(gene_info_list)

# Specify the correct file path for the Excel file
output_file = r'C:/Users/GenepoweRx_Madhu/Downloads/biotype_final_csq_specific_gene.xlsx'
gene_info.to_excel(output_file, sheet_name='Genes', index=False)

print(f"Genes present in VCF file have been saved to {output_file}")
