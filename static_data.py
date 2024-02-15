from app.business_logic.utilities.db_config import DevConfig
from numpy import average

CLN_SIG_PROFILES = ["Hormonal Health Profile", "GastroIntestinal Profile", "Organ Health Profile",
                    "Nutrigenomics Profile", "Fitness Profile", "Blood Factors", "Neuromuscular Factors"]
INSERT_TASK_ARG = "insert_data"
RETRIEVE_TASK_ARG = "get_info"
ALL_INFO_ARG = "info_all"
PROFILE_INFO_ARG = "info_profile"
CATEGORY_INFO_ARG = "info_category"
INFO_FOR_BASED_DATA = "info_for_data"
INFO_FOR_CAT_AND_INFO = "info_for_fgx"
INFO_VIT_MIN = "info_vitamins_minarals"
INSERT_IGX = "insert_igx"
ONLY_CATEGORY_INFO_ARG = "search_by_category"
DISTINCT_PROFILE_ARG = "distinct_profiles"
DISTINCT_RSIDS_ARG = "distinct_rsids"
DISTINCT_GENE_ARG = "distinct_genes"
UPDATE_SYSTEMS_SCORES = "update_systems_scores"
UPDATE_PROFILE = "update_profile"
FITNESS_DATA = "fitness_info"
DATA_BY_IDS = "data_by_ids"
DATA_BY_GENES = "data_by_genes"
ADD_NEW_SCORES = "add_new_scores"
NUTRI_DATA = "nutri_data"
UPDATE_SCORES = "update_scores"
SEARCH_BY_RSIDS = "search_by_rsids"
DISTINCT_CAT_WITH_PROFILE_ARG = "distinct_cat_with_profile"
TRANSCELL_TASK_ARG = "transcell"
WORD_REPORT = "report"
REGULAR_TASK_ARG = "regular"
ANALYTICS_TASK_ARG = "analytics"
SECTIONS_PROFILES_ARG = "sections"
SAMPLE_INFO_AIGPRX = "sample_info"
GENE_INFO_ARG = "gene_info"
VIT_MIN_BY_PROFILE = "by_profile"
DATA_BY_DISEASE = "by_disease"
MEDICAL_DATA = "medical_data"
NGX_DATA = "ngx_data"
FGX_DATA = "fgx_data"
IGX_DATA = "igx_data"
MEDICAL_CONCERN = "medical_concern"
ONE_TIME_UPDATE = "one_time_update"

# These are the profile names present in database.
# If database changes please check this for any update, addition or deletion.

CLINICAL_DISEASE_NAME_COLUMN = "CLNDN"
CONSEQUENCE_COLUMN = "Consequence"
REVIEW_STATUS_COLUMN = "CLNREVSTAT"
ZYGOSITY_COLUMN = "Zygosity"
GENE_NAME_COLUMN = "Gene Name"
RSID_COLUMN = "Rsid"
CASE_SENSITIVE_ATTRIBUTE = "$cs"

# Keyword Lib static content - change here for any updates
KL_KEYWORD_COLUMN = "keyword"
KL_PROFILE_COLUMN = "profile"
KL_CATEGORY_COLUMN = "category"

# FILEs and DIRs
conf = DevConfig().get_config()
BASE_DIR = conf.get('BASE_DIR')
VCF_DIR = f"{BASE_DIR}/vcf"
CLN_SIG_DB_DIR = f"{BASE_DIR}/cln_sig_db"
CLINVAR_DB_DIR = f"{BASE_DIR}/clinvar_db"
OUTPUT_DIR = f"{BASE_DIR}/clinical_significance"
GENE_SCORE_DIR = f"{CLN_SIG_DB_DIR}/Genes_scoring_FIles"
SCORING_CHART_DIR = f"{BASE_DIR}/scoring_chart"
WORD_REPORT_DIR = f"{BASE_DIR}/word_reports"
IMAGES_DIR = f"{GENE_SCORE_DIR}/images"
# TODO Need to change this back to previous one once the run is over.
#OUTPUT_DIR = f"{conf.get('BASE_DIR')}/for_research"
CLN_SIG_DB_FILE = "Gene_library_parsable.xlsx"
DUMMY_FILE_FOR_REPORT = "dummy_for_report.xlsx"
CLINVAR_DB_FILE = "clinvar.vcf"
INDEL_OUTPUT_SUFFIX = "_indel_clinical_significance.xlsx"
VCF_PROCESSED_FILE = "_vcf_processed.xlsx"
INDEL_SHEET_NAME = "Cln Sig Indel"
SNP_VCF_SUFFIX = "_final.vcf"
INDEL_VCF_SUFFIX = "_annotated_indel.vcf"
INHERITANCE_FILE = "clinvar_dbsnp_retireved_data_Genomic data file.xlsx"
SCORING_FILE_SUFFIX = "_Scoring_chart.xlsx"
WORD_FILE_SUFFIX = "_report.docx"

CLASSIFICATION = ['Pathogenic', 'Likely_pathogenic', 'Benign', 'Likely_benign', 'Uncertain_significance',
                  'drug_response', 'association', 'risk_factor', 'protective', 'Affects']

REGULAR_OUTPUT_HEADER = ['Gene Name', 'Rsid', "Final Score", 'Consequence', 'Variant Consequence Score', 'Zygosity',
                         'Zygosity Score',
                         'CLNDN', 'Clinical Consequence Score', 'Weightage', 'Clinical System', 'Clinical Consequence',
                         'Clinical Significance', 'Clinical Significance Score', 'CLNREVSTAT', 'Review Status Score',
                         'Clinvar ID', 'Biotype', 'Strand',
                         'Protein Position and Amino Acid', 'Codons', 'Chromosome', 'Exon', 'Position', 'Allele',
                         'Impact', "Impact Score", 'Existing  Variations', "final_scores_split", "inheritance pattern",
                         "word doc section", "group", "One Liner for word doc"]
CLUBBED_OUTPUT_HEADER = ['Gene Name - Rsid', "Final Score", 'Consequence', 'Variant Consequence Score', 'Zygosity',
                         'Zygosity Score',
                         'CLNDN', 'Clinical Consequence Score', 'Weightage', 'Clinical System', 'Clinical Consequence',
                         'Clinical Significance', 'Clinical Significance Score', 'CLNREVSTAT', 'Review Status Score',
                         'Clinvar ID', 'Biotype', 'Strand',
                         'Protein Position and Amino Acid', 'Codons', 'Chromosome', 'Exon', 'Position', 'Allele',
                         'Impact', "Impact Score", 'Existing  Variations', "final_scores_split", "inheritance pattern",
                         "word doc section", "group", "One Liner for word doc"]

TRANSCELL_OUTPUT_HEADER = ['Gene Name', 'Rsid', "Final Score", 'Consequence', 'Variant Consequence Score', 'Zygosity',
                           'Zygosity Score',
                           'CLNDN', 'Clinical System',
                           'Clinical Significance', 'Clinical Significance Score', 'CLNREVSTAT', 'Review Status Score',
                           'Clinvar ID', 'Biotype', 'Strand',
                           'Protein Position and Amino Acid', 'Codons', 'Chromosome', 'Exon', 'Position', 'Allele',
                           'Impact', "Impact Score", 'Existing  Variations', "final_scores_split",
                           "inheritance pattern",
                           "word doc section", "group", "One Liner for word doc"]

CLUBBED_HEADER_GENE_RSID = "Gene Name - Rsid"

# Hormonal Health, Gastro and organ health need to be generated in both the ways systems as well as nutrigenomics.

SYSTEM_PROFILES_CATEGORIES = {
    "Hormonal Health Profile": ["Hypothalamus & Pituitary", "Pineal gland", "Thyroid & Parathyroid",
                                "Adrenal gland", "Ovaries", "Testes", "Kidneys", "Pancreas"],
    "GastroIntestinal Profile": ["Esophageal Disorders", "Stomach/Gastric Disorders",
                                 "Intestinal Disorders/ Colorectal Disorders", "Pancreatic Disorders",
                                 "Biliary Disorders"],
    "Organ Health Profile": ["Diabetes Mellitus", "Hypertension and Cardiovascular health", "Kidney health",
                             "Lung health", "Bone health", "Skin / Derm health", "Eye / Ear health",
                             "Parkinsons", "Alzheimers", "Anxiety / Depression/ Schizophrenia",
                             "Cancer predisposition", "Mitochondrial disorders", "Reproductive health -Women",
                             "Reproductive health -Men", "MISCELLANEOUS"]
}

HORMONAL_HEALTH_KEYWORD_CATEGORIES = ["Hypothalamus & Pituitary", "Pineal gland", "Thyroid & Parathyroid",
                                      "Adrenal gland", "Ovaries", "Testes", "Kidneys", "Pancreas"]

GASTROINTESTINAL_KEYWORD_CATEGORIES = ["Esophageal Disorders", "Stomach/Gastric Disorders",
                                       "Intestinal Disorders/ Colorectal Disorders", "Pancreatic Disorders",
                                       "Biliary Disorders"]

ORGAN_HEALTH_KEYWORD_CATEGORIES = ["Diabetes Mellitus", "Hypertension and Cardiovascular health", "Kidney health",
                                   "Lung health", "Bone health", "Skin / Derm health", "Eye / Ear health",
                                   "Parkinsons", "Alzheimers", "Anxiety / Depression/ Schizophrenia",
                                   "Cancer predisposition", "Mitochondrial disorders", "Reproductive health -Women",
                                   "Reproductive health -Men", "MISCELLANEOUS"]

NUTRIGENOMICS_CATEGORIES = ["Carbohydrate Metabolism", "Lipid Metabolism", "Protein Metabolism",
                            "Vitamin and Minerals Disorders", "Other Metabolic Disorders"]

NEUROMUSCULAR_CATEGORIES = ["Spinal muscular atrophies", "Peripheral nerve disorders", "Neuromuscular junction",
                            "Metabolic diseases of the muscle", "Other rare myopathies"]

FITNESS_CATEGORIES = ["Muscle /Joint / Connective Tissue Factors", "IGX Factors"]

BLOOD_FACTOR_CATEGORIES = ["Blood Factor"]

OMICS_PROFILES = ["Hormonal Health Profile", "GastroIntestinal Profile", "Nutrigenomics Profile", "Fitness Profile",
                  "Blood Factors", "Neuromuscular Factors"]
SYSTEM_PROFILE_FILE_NAME_PART = "Systems"
# Word Documents templates info
PREGEN = "pregen"
SCELL = "scell"

PATHOGENIC_KEYS = ["patho_gene", "mutations", "score", "phy_reco"]

# MB scores
# Variant type: VT
VT_SNV = 50
VT_INDEL = 50
VT_CNV = 75

# clinical significnace type: CST
CST_PATHOGENIC = 75
CST_RISK = 50
CST_PROTECTIVE = 15
CST_BENIGN = 25
CST_AFFECTS = 10
CST_ASSOCIATIONS = 10
CST_UNCERTAIN = 10

# zygosity : Z
Z_HOMO = 100
Z_HETERO = 50

# consequence: CON
CON_SYNONYMOUS = 10
CON_MISSENSE = 60
CON_FRAMESHIFT = 50
CON_SPLICE = 50
CON_UTR = 25
CON_STOP_LOSS = 50
CON_STOP_GAIN = 75
CON_NONSENSE = 75
CON_INFRAME_DELETION = 25

# Review status: RS
RS_PRACTICE = 200
RS_EXPERT = 75
RS_MULTI_NO_CONFLICT = 25

# modifier effect: MODI
MODI_LOW = 25
MODI_MODERATE = 50
MODI_HIGH = 100

MB_LOW_HIGH = MODI_HIGH + CON_STOP_GAIN + Z_HOMO + CST_PATHOGENIC + VT_SNV
HIGHEST_MB_SCORE_POSSIBLE = round(float((4 * MB_LOW_HIGH + 6 * RS_PRACTICE) / 10))

MB_LOW_LOW = MODI_LOW + Z_HETERO + CON_SYNONYMOUS + 0 + VT_SNV
LOWEST_MB_SCORE_POSSIBLE = round(float((4 * MB_LOW_LOW + 6 * RS_MULTI_NO_CONFLICT) / 10))

# Normalization Formula: 0 + float(x-LOWEST_MB_SCORE_POSSIBLE)* float(10 - 0) / (500 - 160)

# polyphen or sift scores need to be added next.

# consequence being considered
CONSEQUENCE_BEING_CONSIDERED = ["missense variant", "frameshift variant", "inframe deletion",
                                "inframe insertion", "inframe indel", "stop lost", "stop gain",
                                "nonsense", "splice donor variant", "splice acceptor variant"]

INTRON_CONSEQUENCES = ["intron variant", "non-coding transcript variant", "5 prime UTR variant",
                       "3 prime UTR variant", "initiatior codon variant",
                       "genic upstream transcript variant", "genic downstream transcript variant",
                       "no sequence alteration", "synonymous variant"]

KEYWORDS_INTRON = ["intron", "transcript", "UTR", "codon", "upstream", "downstream", "sequence", "synonymous"]
KEYWORDS_EXON = ["missense", "frame", "stop", "nonsense", "splice", "indel"]

# Weights for weighted average of clinical score and MB score
WEIGHTS = [30, 70]
REVERSE_WEIGHTS = [70, 30]
LOWEST_WEIGHTED_AVERAGE = round(average([LOWEST_MB_SCORE_POSSIBLE, 0], weights=WEIGHTS), 2)
HIGHEST_WEIGHTED_AVERAGE = round(average([HIGHEST_MB_SCORE_POSSIBLE, 10], weights=WEIGHTS), 2)

LOWEST_REVERSE_WEIGHTED_AVERAGE = round(average([LOWEST_MB_SCORE_POSSIBLE, 0], weights=REVERSE_WEIGHTS), 2)
HIGHEST_REVERSE_WEIGHTED_AVERAGE = round(average([HIGHEST_WEIGHTED_AVERAGE, 10], weights=REVERSE_WEIGHTS), 2)

SYSTEM_CLASSIFICATION = ["Cardiovascular Health", "ENT", "Gastrointestinal Health", "Hormonal Health", "Lungs Health",
                         "Rare Disorder", "Skin Health", "blood factors", "cancers", "cholesterol disorders",
                         "diabetes", "eye health", "fitness genomics", "immunity", "metabolic health",
                         "musculoskeletal system", "nervous system", "nutrigenomics", "obesity", "renal health"]
INHERITANCE_COLUMNS = ["Genes", "RSID", "CLNDN (Main) - No modifications Done", "CLNDN",
                       "One line CLNDN for dummy data",
                       "Section Name", "Inheritance Pattern_New", "Group (Adult/Pediatric/Neonatal/Both)",
                       "CLNSIG", "CLNREVSTAT", "Molecular_consequences", "Clinvar_ID", "Chromosome",
                       "Position", "Alleles"]
INHERITANCE_DB_COLUMNS = ["genes", "rsid", "clndn", "one_liner", "section_name", "inheritance_pattern",
                          "grouped_under", "clnsig", "review_status", "molecular_consequence", "clinvar_id",
                          "chromosome", "position_of_mutant", "alleles"]

SCELLCARE_CANCER_GENES = ["MTOR", "MTHFR", "DPYD", "MSH2", "MSH6", "FBXO11", "MLH1",
                          "PIK3CA", "PMS2", "AIMP2", "EGFR", "ABCB1", "BRCA2", "BRCA1",
                          "MYH7", "MAP2K1", "CDH1", "NQO1"]

BOTH_HOM_HET_INHERITANCE = ["AD", "AD, AR", "XLD, XLR", "AD, DD", "AD, AR, DD"]
GENDER_BASED_INHERITANCE_MALE = {
    "XLD": "Heterozygous",
    "XLR": "Heterozygous",
    "XL": "Heterozygous"
}
GENDER_BASED_INHERITANCE_FEMALE = {
    "XLD": "both",
    "XLR": "Homozygous",
    "XL": "both"
}

ONLY_HOM_INHERITANCE = ["AR", "-", "None", "AR, DD", "AR, DR"]

# Clinical Genomics constants
CLINICAL_GENOMICS_DIR = f"{BASE_DIR}/cln_sig_db/clinical_genomics"
UNIQUE_CLNDN_FILE = f"{CLINICAL_GENOMICS_DIR}/K&H genomic db_Clinical Conditions.xlsx"
MAIN_CLINICAL_GENOMICS_FILE = f"{CLINICAL_GENOMICS_DIR}/inheritance_clinvar_transcell.xlsx"
INHERITANCE_FILES_PREFIX = "inheritance_clinvar_transcell_part"

CONDITION_LOW_TO_MILD = "Low to Mild"
CONDITION_MILD = "Mild"
CONDITION_MILD_TO_MODERATE = "Mild to Moderate"
CONDITION_MODERATE = "Moderate"
CONDITION_MODERATE_TO_HIGH = "Moderate to High"
CONDITION_HIGH = "High"
CONDITION_LOW = "Low"
CONDITION_NIL = "no mutations"
CONDITION_IS_CONCERN = "is_concern"
FATIGUE_CONDITION_INFO = "fatigue_info"
FOOD_INFO = "food_info"
IMMUNITY_CAT = "Immunity"
NON_MUTUALLY_EXCLUSIVE = [CONDITION_IS_CONCERN, FATIGUE_CONDITION_INFO, FOOD_INFO]

MEDICAL_CONDITION_MAP = {
    1: "Diabetes",
    2: "High Blood pressure",  # High Blood pressure and Low Blood pressure come one below other.
    3: "Coronary Artery Disease",
    4: "Arrhythmia",
    5: "Heart Failure- Dilated Cardiomyopathy, Restrictive Cardiomyopathy",
    6: "Cholesterol disorders",
    7: "Hypertriglyceridemia",
    8: "Thyroid Disorders- Hypothyroidism, Hyperthyroidism",
    9: "Anemia- Microcytic, Hemolytic",
    10: "Predisposition to Blood clots- Thrombophilia",
    11: "Bleeding Disorders",
    12: "Parkinson’s Disease",
    13: "Alzheimer’s Disease",
    14: "Migraines, Headaches",
    15: "Seizures",
    16: "Inflammatory bowel disease- Crohn’s, Ulcerative colitis",
    17: "Respiratory Allergies",
    18: "Food Allergies",
    19: "Liver Disorders",
    20: "Gall bladder disorders",
    21: "Pancreatic Disorders",
    22: "Nephrotic Syndrome (Focal Segmental Glomerulosclerosis, Membranous nephropathy, Minimal Change Disease)",
    23: "Interstitial Nephritis, Tubulo interstitial Disease",
    24: "Renal Stones- Calcium Oxalate stones, Cystine stones, Uric Acid Stones",
    25: "Dry skin, eczema",
    26: "Skin Allergies",
    27: "Vitiligo",
    28: "Osteoporosis",
    29: "Degenerative Joint Disease, Cartilage degeneration",
    30: "Muscular dystrophy, atrophy",
    31: "Fatigue",  # Need to consider 9, 28, 29, 30 for fatigue as well. ( put in the logic)
    32: "Mood Disorders- Anxiety, Schizophrenia, Depression",
    33: "Urticaria", # Condition Specific
    34: "Essential tremors", # Condition Specific
    35: "Renal Disorders",
    36: "Sinusitis, Dust Allergy (Ciliary dykinesia, Hyper IgE syndrome, Angioedma, Chroinc granulomatous)",
    37: "Obesity",
    38: "Skin Health",
    39: "Eye and Ear health",
    40: "Low Blood Pressure",
    41: "Gastritis", # Same as IBD
    42: "Glaucoma"  # Still do not have any info yet.
}
