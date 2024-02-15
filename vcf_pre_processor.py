import numpy as np
import sys

from app.business_logic.utilities.file_reader import vcf_reader
from collections import OrderedDict
from app.business_logic.utilities import static_data as sd
from app.business_logic.utilities import output_formatter as ofmt
from app.business_logic.utilities import file_writer as wt
import pandas as pd



class VcfPreProcess:
    def __init__(self, sample, snp_vcf, indel_vcf):
        self.sample = sample
        self.snp_vcf = snp_vcf
        self.indel_vcf = indel_vcf

    def start_process(self):
        # indel_records = vcf_reader(self.indel_vcf)
        # pre-process Indel records

        # TODO Use khvcf sent in mail to parse the vcf and get the snp records.
        snp_records = vcf_reader(self.snp_vcf)
        # pre-process snp records
        list_of_rows = preprocess_snp_records(snp_records)
        wt.write_pre_processed_vcf(list_of_rows, self.sample)
        #haplotype_rsid_zygosity(snp_records, self.sample)


def haplotype_rsid_zygosity(vcf_records, sample):
    rsid_data = pd.read_excel("C:/Users/win/code/data/cln_sig_db/halplotype_info_sheet.xlsx")

    rsid_data = rsid_data[rsid_data.rsid != "-"].drop_duplicates()
    rsid_list = rsid_data.rsid.tolist()
    print(rsid_list)
    print(rsid_data)
    final_data = []
    vcf_info = []
    header = vcf_records.infos['CSQ'].desc.split(':')[1].split("|")
    for record in vcf_records:
        try:
            rsid = record.ID.strip()
        except Exception:
            rsid = "-"
        # Zygocity
        if record.INFO['HET'] == 1:
            zygocity = 'Heterozygous'
        elif record.INFO['HOM'] == 1:
            zygocity = 'Homozygous'
        else:
            zygocity = 'Unknown'
        if rsid in rsid_list:
            rsid_data.loc[rsid_data.rsid == rsid, "Zygosity"] = zygocity

    wild = rsid_data[(rsid_data.rsid.notnull()) & (rsid_data.Zygosity.isnull())]
    unwild = rsid_data[(rsid_data.rsid.notnull()) & (rsid_data.Zygosity.notnull())]
    wild["Zygosity"] = "Wild Type"
    print(wild)
    print(unwild)
    frames = [unwild, wild]
    final_info = pd.concat(frames, ignore_index=True, sort=False)

    final_info.to_excel(f"C:/Users/win/code/data/pharma_db/haplotype_rsid_zygosity_output_{sample}.xlsx")


def create_dummy_row(final_data, rsid_list, gene_rsid):
    wild_type = []
    fd_rsids = [r["Rsid"] for r in final_data]
    non_selected_rsids = [r for r in rsid_list if r not in fd_rsids]
    for rsid in non_selected_rsids:
        row = OrderedDict({
            sd.REGULAR_OUTPUT_HEADER[0]: get_gene(rsid, gene_rsid),
            sd.REGULAR_OUTPUT_HEADER[1]: rsid,
            sd.REGULAR_OUTPUT_HEADER[2]: "-",
            sd.REGULAR_OUTPUT_HEADER[3]: 0,
            sd.REGULAR_OUTPUT_HEADER[4]: "Wild Type",
            sd.REGULAR_OUTPUT_HEADER[5]: 0,
            sd.REGULAR_OUTPUT_HEADER[6]: "-",
            sd.REGULAR_OUTPUT_HEADER[7]: 0,
            sd.REGULAR_OUTPUT_HEADER[8]: "NA",
            sd.REGULAR_OUTPUT_HEADER[9]: "NA",
            sd.REGULAR_OUTPUT_HEADER[10]: "-",
            sd.REGULAR_OUTPUT_HEADER[11]: 0,
            sd.REGULAR_OUTPUT_HEADER[12]: "-",
            sd.REGULAR_OUTPUT_HEADER[13]: 0,
            sd.REGULAR_OUTPUT_HEADER[14]: "-",
            sd.REGULAR_OUTPUT_HEADER[15]: "-",
            sd.REGULAR_OUTPUT_HEADER[16]: "-",
            sd.REGULAR_OUTPUT_HEADER[17]: "-",
            sd.REGULAR_OUTPUT_HEADER[18]: "-",
            sd.REGULAR_OUTPUT_HEADER[19]: "-",
            sd.REGULAR_OUTPUT_HEADER[20]: "-",
            sd.REGULAR_OUTPUT_HEADER[21]: "-",
            sd.REGULAR_OUTPUT_HEADER[22]: "-",
            sd.REGULAR_OUTPUT_HEADER[23]: "-",
            sd.REGULAR_OUTPUT_HEADER[24]: "-"
        })
        wild_type.append(row)
    return wild_type

def hap_info(vcf_info):
    hap_data = []
    # Read haplotype_info_sheet

    # separate out the rsids into a list

    # for each record in vcf list check if its in the list and add zygosity to hapdata


def get_gene(rsid, rsid_df):
    gene = rsid_df[rsid_df.rsid ==rsid]
    print(gene)
    return gene['Gene'].values[0]


def preprocess_snp_records(vcf_records):
    vcf_rows = []
    if vcf_records:
        print(type(vcf_records))
    header = vcf_records.infos['CSQ'].desc.split(':')[1].split("|")
    for record in vcf_records:
        list_of_rows = []
        chromosome = record.CHROM
        position = record.POS
        rsid = record.ID
        allele = str(record.alleles)
        allele = allele.replace("'", " ")
        # Zygocity
        if record.INFO['HET'] == 1:
            zygocity = 'Heterozygous'
        elif record.INFO['HOM'] == 1:
            zygocity = 'Homozygous'
        else:
            zygocity = 'Unknown'
        try:
            gene_name = record.INFO['GENEINFO'].split(':')[0]
        except KeyError:
            # print(f"Gene Information not present for this record {record}")
            gene_name = 'N/A'
        for csq_parts in record.INFO['CSQ']:
            temp_dict = {header[i]: csq_parts.split('|')[i] for i in range(len(header))}
            if temp_dict["Amino_acids"] and temp_dict["Protein_position"] \
                    and len(temp_dict["Amino_acids"].split('/')) > 1:
                protein_pos_amino_acid = f'{temp_dict["Amino_acids"].split("/")[0]}' \
                                         f'{temp_dict["Protein_position"]}' \
                                         f'{temp_dict["Amino_acids"].split("/")[1]}'
            else:
                protein_pos_amino_acid = 'N/A'
            row = OrderedDict({
                sd.REGULAR_OUTPUT_HEADER[0]: gene_name,
                sd.REGULAR_OUTPUT_HEADER[1]: rsid,
                sd.REGULAR_OUTPUT_HEADER[2]: ofmt.format_consequence_disease(temp_dict['Consequence']),
                sd.REGULAR_OUTPUT_HEADER[3]: 0,
                sd.REGULAR_OUTPUT_HEADER[4]: zygocity,
                sd.REGULAR_OUTPUT_HEADER[5]: 0,
                sd.REGULAR_OUTPUT_HEADER[6]: ofmt.format_consequence_disease(temp_dict['ClinVar_CLNDN']),
                sd.REGULAR_OUTPUT_HEADER[7]: 0,
                sd.REGULAR_OUTPUT_HEADER[8]: "NA",
                sd.REGULAR_OUTPUT_HEADER[9]: "NA",
                sd.REGULAR_OUTPUT_HEADER[10]: ofmt.format_review_status_cln_sig(temp_dict['ClinVar_CLNSIG']),
                sd.REGULAR_OUTPUT_HEADER[11]: 0,
                sd.REGULAR_OUTPUT_HEADER[12]: ofmt.format_review_status_cln_sig(temp_dict['ClinVar_CLNREVSTAT']),
                sd.REGULAR_OUTPUT_HEADER[13]: 0,
                sd.REGULAR_OUTPUT_HEADER[14]: temp_dict['ClinVar'],
                sd.REGULAR_OUTPUT_HEADER[15]: temp_dict['BIOTYPE'],
                sd.REGULAR_OUTPUT_HEADER[16]: temp_dict['STRAND'],
                sd.REGULAR_OUTPUT_HEADER[17]: protein_pos_amino_acid,
                sd.REGULAR_OUTPUT_HEADER[18]: temp_dict['Codons'],
                sd.REGULAR_OUTPUT_HEADER[19]: chromosome,
                sd.REGULAR_OUTPUT_HEADER[20]: temp_dict['EXON'],
                sd.REGULAR_OUTPUT_HEADER[21]: position,
                sd.REGULAR_OUTPUT_HEADER[22]: allele,
                sd.REGULAR_OUTPUT_HEADER[23]: temp_dict['IMPACT'],
                sd.REGULAR_OUTPUT_HEADER[24]: temp_dict['Existing_variation']
            })
            list_of_rows.append(row)
        #list_of_rows = [i for n, i in enumerate(list_of_rows) if i not in list_of_rows[n + 1:]]
        vcf_rows.extend(list_of_rows)
    return vcf_rows


if __name__ == "__main__":
    print("This is running as standalone")
    sample = "KHHSPTGPPGX6"
    snp_vcf = f"{sd.VCF_DIR}/{sample}{sd.SNP_VCF_SUFFIX}"
    indel_vcf = f"{sd.VCF_DIR}/{sample}{sd.INDEL_VCF_SUFFIX}"
    ppv = VcfPreProcess(sample, snp_vcf, indel_vcf)
    ppv.start_process()
