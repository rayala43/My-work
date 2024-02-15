from app.business_logic.utilities import static_data as sd
from collections import OrderedDict, defaultdict
import re
from numpy import average


# String formatting of values from vcf.
# TODO need to introduce flask logging
def format_consequence_disease(clndn_str):
    input_string = clndn_str.replace('_', ' ')
    input_string = format_not_specified(input_string)
    input_string = format_not_provided(input_string)
    input_string = format_none_provided(input_string)
    input_string = format_to_compliance(input_string, "resistance to")
    input_string = format_to_compliance(input_string, "susceptibility to")
    input_string = format_to_compliance(input_string, "protective against")
    input_string = input_string.replace('_', ' ').replace('&', ', ').replace('|', ', ').strip().rstrip(",")

    return input_string


def format_review_status_cln_sig(revstat):
    return revstat.replace('_', ' ').replace('&', ', ')


def format_none_provided(input_string):
    if input_string.replace("none provided", "") \
            and input_string.replace("none provided", "") != "" \
            and input_string.replace("None provided", "") and input_string.replace("None provided", "") != "":
        input_string = input_string.replace("&none provided", '').replace("none provided&", '')
        input_string = input_string.replace("&None provided", "").replace("None provided&", "")
        input_string = input_string.replace("|none provided", '').replace("none provided|", '')
        input_string = input_string.replace("|None provided", "").replace("None provided|", "")
        input_string = re.sub(r"&None specified$", "", input_string, flags=re.IGNORECASE)


    return input_string


def format_not_specified(input_string):
    if input_string.replace("not specified", "") \
            and input_string.replace("not specified", "") != "" \
            and input_string.replace("Not specified", "") and input_string.replace("Not specified", "") != "":
        input_string = input_string.replace("&not specified", '').replace("not specified&", '')
        input_string = input_string.replace("&Not specified", "").replace("Not specified&", "")
        input_string = input_string.replace("|Not specified", "").replace("Not specified|", "")
        input_string = input_string.replace("|not specified", "").replace("not specified|", "")
        input_string = re.sub(r"&Not specified$", "", input_string, flags=re.IGNORECASE)

    return input_string


def format_not_provided(input_string):
    if input_string.replace("not provided", "") \
            and input_string.replace("not provided", "") != "" \
            and input_string.replace("Not provided", "") and input_string.replace("Not provided", "") != "":
        input_string = input_string.replace("&not provided", "").replace("not provided&", "")
        input_string = input_string.replace("&Not provided", "").replace("Not provided&", "")
        input_string = input_string.replace("|not provided", "").replace("not provided|", "")
        input_string = input_string.replace("|Not provided", "").replace("Not provided|", "")
        input_string = re.sub(r"&Not provided$", "", input_string, flags=re.IGNORECASE)

    return input_string


def format_to_compliance(input_string, search_string):
    try:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        s_string = re.search(search_string, input_string, flags=re.IGNORECASE).group()

        indices = []
        i = 0
        for p in parts:
            if p.lower() == s_string.lower():
                indices.append(i)
            i += 1
        for index_of_susceptibility in indices:
            if index_of_susceptibility == 0:
                input_string = input_string.replace(f"{s_string}&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
        input_string = ", ".join(parts).replace(f"{s_string},", f"{s_string}")
        return input_string
    except Exception:
        pass
    return input_string


def format_susceptibility_to(input_string):
    if "susceptibility to" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('susceptibility to')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("susceptibility to&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("susceptibility to,", "susceptibility to")
        except Exception as es:
            # print(es)
            pass
    elif "Susceptibility to" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('Susceptibility to')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("Susceptibility to&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("Susceptibility to,", "Susceptibility to")
        except Exception as es:
            # print(es)
            pass
    elif "SUSCEPTIBILITY TO" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('SUSCEPTIBILITY TO')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("SUSCEPTIBILITY TO&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("SUSCEPTIBILITY TO,", "SUSCEPTIBILITY TO")
        except Exception as es:
            # print(es)
            pass
    return input_string


def format_resistance_to(input_string):
    if "resistance to" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('resistance to')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("resistance to&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("resistance to,", "resistance to")
        except Exception as es:
            # print(es)
            pass
    elif "Resistance to" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('Resistance to')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("Resistance to&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("Resistance to,", "Resistance to")
        except Exception as es:
            # print(es)
            pass
    elif "RESISTANCE TO" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('RESISTANCE TO')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("RESISTANCE TO&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("RESISTANCE TO,", "RESISTANCE TO")
        except Exception as es:
            # print(es)
            pass
    return input_string


def format_protective_against(input_string):
    if "protective against" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('protective against')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("protective against&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("protective against,", "protective against")
        except Exception as es:
            # print(es)
            pass
    elif "Protective against" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('Protective against')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("Protective against&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("Protective against,", "Protective against")
        except Exception as es:
            # print(es)
            pass
    elif "PROTECTIVE AGAINST" in input_string:
        parts = input_string.split('&')
        parts = [x.strip() for x in parts]
        try:
            index_of_susceptibility = parts.index('PROTECTIVE AGAINST')
            if index_of_susceptibility == 0:
                input_string = input_string.replace("PROTECTIVE AGAINST&", "")
            else:
                parts[index_of_susceptibility - 1], parts[index_of_susceptibility] = parts[index_of_susceptibility], \
                                                                                     parts[index_of_susceptibility - 1]
                input_string = ", ".join(parts).replace("PROTECTIVE AGAINST,", "PROTECTIVE AGAINST")
        except Exception as es:
            # print(es)
            pass
    return input_string


# Formatting the rows (pooling them based on rsid) to confirm to output excel required.
def pool_with_rsids(rows, is_clubbed, cln_score_required):
    pooled_rows = []
    cond_dict = defaultdict(list)
    for row in rows:
        cond_dict[row[sd.CLINICAL_DISEASE_NAME_COLUMN]].append(row)
    gene_dicts = []
    for k, v in cond_dict.items():
        gd = defaultdict(list)
        for s in v:
            gd[s[sd.GENE_NAME_COLUMN]].append(s)
        gene_dicts.append(gd)
    gene_dicts = list(gene_dicts)
    zyg_dicts = []
    for g in gene_dicts:
        for k, v in g.items():
            temp_dict = {}
            hom = [x for x in v if x[sd.ZYGOSITY_COLUMN] == "Homozygous"]
            het = [x for x in v if x[sd.ZYGOSITY_COLUMN] == "Heterozygous"]
            temp_dict[f"{k}_hom"] = hom
            temp_dict[f"{k}_het"] = het
            zyg_dicts.append(temp_dict)
    for r in zyg_dicts:
        # For dict get the values and club them to a single row
        for k, v in r.items():
            if v and v != []:
                vs = v[0]
                rsids = [rs[sd.RSID_COLUMN] if rs[sd.RSID_COLUMN] is not None else "No_rsid" for rs in v]

                scores_dict = add_scores(v, cln_score_required)
                temp_dict = OrderedDict()
                if is_clubbed:
                    temp_dict[sd.CLUBBED_HEADER_GENE_RSID] = club_gene_rsids(vs[sd.REGULAR_OUTPUT_HEADER[0]], rsids)
                else:
                    temp_dict[sd.REGULAR_OUTPUT_HEADER[0]] = vs[sd.REGULAR_OUTPUT_HEADER[0]]
                    temp_dict[sd.REGULAR_OUTPUT_HEADER[1]] = without_clubbing(rsids)
                temp_dict[sd.REGULAR_OUTPUT_HEADER[2]] = scores_dict["total_final"]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[3]] = vs[sd.REGULAR_OUTPUT_HEADER[3]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[4]] = scores_dict["total_consequence"]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[5]] = vs[sd.REGULAR_OUTPUT_HEADER[5]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[6]] = scores_dict["total_zygosity"]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[7]] = vs[sd.REGULAR_OUTPUT_HEADER[7]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[8]] = vs[sd.REGULAR_OUTPUT_HEADER[8]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[9]] = vs[sd.REGULAR_OUTPUT_HEADER[9]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[10]] = vs[sd.REGULAR_OUTPUT_HEADER[10]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[11]] = vs[sd.REGULAR_OUTPUT_HEADER[11]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[12]] = vs[sd.REGULAR_OUTPUT_HEADER[12]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[13]] = scores_dict["total_cln_sig_score"]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[14]] = vs[sd.REGULAR_OUTPUT_HEADER[14]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[15]] = scores_dict["total_revstat"]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[16]] = vs[sd.REGULAR_OUTPUT_HEADER[16]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[17]] = vs[sd.REGULAR_OUTPUT_HEADER[17]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[18]] = vs[sd.REGULAR_OUTPUT_HEADER[18]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[19]] = vs[sd.REGULAR_OUTPUT_HEADER[19]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[20]] = vs[sd.REGULAR_OUTPUT_HEADER[20]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[21]] = vs[sd.REGULAR_OUTPUT_HEADER[21]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[22]] = vs[sd.REGULAR_OUTPUT_HEADER[22]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[23]] = vs[sd.REGULAR_OUTPUT_HEADER[23]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[24]] = vs[sd.REGULAR_OUTPUT_HEADER[24]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[25]] = vs[sd.REGULAR_OUTPUT_HEADER[25]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[26]] = scores_dict["total_modi"]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[27]] = vs[sd.REGULAR_OUTPUT_HEADER[27]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[28]] = scores_dict["final_split"]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[29]] = vs[sd.REGULAR_OUTPUT_HEADER[29]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[30]] = vs[sd.REGULAR_OUTPUT_HEADER[30]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[31]] = vs[sd.REGULAR_OUTPUT_HEADER[31]]
                temp_dict[sd.REGULAR_OUTPUT_HEADER[32]] = vs[sd.REGULAR_OUTPUT_HEADER[32]]
                pooled_rows.append(temp_dict)
    return pooled_rows


# Helpers for pooling by rsid function
def club_gene_rsids(gene, rsids):
    rss = "\n".join(rsids)
    gr_set = set(rss.split("\n"))
    gr_string = ""
    for gr in gr_set:
        if gr and gr.strip() != "":
            gr_string = f"{gr_string}\n{gr}"
    gene = f"{gene}{gr_string}"
    return gene


def add_scores(rows, cln_scr_required):
    scores_dict = {
        "total_final": 0,
        "total_zygosity": 0,
        "total_consequence": 0,
        "total_modi": 0,
        "total_clinical_score": 0,
        "total_cln_sig_score": 0,
        "total_revstat": 0,
        "final_split": ""
    }
    clinical_score = 0
    for row in rows:
        if row["Clinical Consequence Score"] > clinical_score:
            clinical_score = row["Clinical Consequence Score"]
        scores_dict["total_zygosity"] += row["Zygosity Score"]
        scores_dict["total_consequence"] += row["Variant Consequence Score"]
        scores_dict["total_modi"] += row["Impact Score"]
        if cln_scr_required:
            scores_dict["total_clinical_score"] += row["Clinical Consequence Score"]
        else:
            scores_dict["total_clinical_score"] = 0
        scores_dict["total_cln_sig_score"] += row["Clinical Significance Score"]
        scores_dict["total_revstat"] += row["Review Status Score"]

        # Calculating final score for the row and finding the total final score.
        mb_low_relavance_score = row["Zygosity Score"] + row["Variant Consequence Score"] + row["Impact Score"] + row[
            "Clinical Significance Score"] + sd.VT_SNV
        total_score_mb = round(float((4 * mb_low_relavance_score + 6 * row["Review Status Score"])/10))
        normalized_mb_score = 0 + float(total_score_mb - sd.LOWEST_MB_SCORE_POSSIBLE
                                        ) * float(10 - 0) / (sd.HIGHEST_MB_SCORE_POSSIBLE - sd.LOWEST_MB_SCORE_POSSIBLE)
        row["Final Score"] = normalized_mb_score
        scores_dict["total_final"] += row["Final Score"]
        scores_dict["final_split"] = f"{scores_dict['final_split']} + {normalized_mb_score}"

    tf_average = float(scores_dict["total_final"])/float(len(rows))
    if cln_scr_required:
        weights = [float(1 - (float(clinical_score) / float(100))), float(clinical_score) / float(100)]
        if 5 <= len(rows) < 11:
            tf_average = tf_average + float(0.5)
        elif len(rows) > 10:
            tf_average = tf_average + float(1)
        scores_dict["total_final"] = round(average([tf_average, clinical_score], weights=weights), 2)
    else:
        scores_dict["total_final"] = round(tf_average, 2)

    return scores_dict


def final_score_calculation_with_decision():
    pass


def without_clubbing(rsids):
    rs_set = set(rsids)
    rs_string = ""
    for rs in rs_set:
        if rs_string and rs_string != "":
            rs_string = f"{rs_string}, {rs}"
        else:
            rs_string = rs
    return rs_string


def consolidate_rsids_conditions(rows):
    final_rows = []
    rsid_rows = {}
    for row in rows:
        if row["Rsid"] and row["Rsid"].strip() and row["Rsid"].strip() != "":
            if not rsid_rows or rsid_rows == {}:
                rsid_rows[row["Rsid"]] = []
                rsid_rows[row["Rsid"]].append(row)
            elif row["Rsid"] in rsid_rows:
                rsid_rows[row["Rsid"]].append(row)
            else:
                rsid_rows[row["Rsid"]] = []
                rsid_rows[row["Rsid"]].append(row)

    for k, v in rsid_rows.items():

        temp_row = {}
        for val_row in v:
            if not temp_row or temp_row == {}:
                temp_row = val_row
            else:
                parts_val = val_row["Consequence"].split(",")
                parts_temp = temp_row["Consequence"].split(",")
                parts_temp.extend(parts_val)
                parts_temp = [x.strip() for x in parts_temp]
                con_set = list(OrderedDict.fromkeys(parts_temp))
                temp_row["Consequence"] = ",".join(con_set)
        ps = temp_row["Consequence"].split(",")
        temp_row["Consequence"] = "\n".join(ps)
        if temp_row["Rsid"] == "rs2919360":
            print("We are Here as well")
        final_rows.append(temp_row)

    return final_rows


# Creating specific rows for formatted excel output
def create_category_row(category):
    row = {}
    for header in sd.CLUBBED_OUTPUT_HEADER:
        if header == sd.CLINICAL_DISEASE_NAME_COLUMN:
            row[header] = category
        else:
            row[header] = " "
    return row


def no_info_row():
    row = OrderedDict()
    for header in sd.CLUBBED_OUTPUT_HEADER:
        if header == sd.CLINICAL_DISEASE_NAME_COLUMN:
            row[header] = "No Homozygous Missense Variant found"
        elif header == sd.CONSEQUENCE_COLUMN:
            row[header] = "No info"
        else:
            row[header] = " "
    return row


# temporary category and no info rows for systems files
def create_category_row_temp(category):
    row = {}
    for header in sd.CLUBBED_OUTPUT_HEADER:
        if header == sd.CLINICAL_DISEASE_NAME_COLUMN:
            row[header] = category
        else:
            row[header] = " "
    return row


def no_info_row_temp():
    row = OrderedDict()
    for header in sd.CLUBBED_OUTPUT_HEADER:
        if header == sd.CLINICAL_DISEASE_NAME_COLUMN:
            row[header] = "No Homozygous Missense Variant found"
        elif header == sd.CONSEQUENCE_COLUMN:
            row[header] = "No info"
        else:
            row[header] = " "
    return row
