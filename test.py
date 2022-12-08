import numpy as np
import pandas as pd

CELLINFO_CCLE_PATH = './raw_data/cellinfo_CCLE.csv'
TCGA_DRUG_TREATMENT_PATH = './raw_data/tcga_drug_treatment.csv'
DRUG_INFO_PATH = './raw_data/drug_info.csv'
BHK_DRUG_INFO_PATH = './raw_data/BHK_drug_info.csv'


def get_dictionary(df, col_name0, col_name1):
    result = df.groupby(col_name0)[col_name1].apply(list).to_dict()
    return result


def search_in_list(_v, _list):
    i = 0
    for row in _list:
        j = 0
        for unit in row:
            if _v in str(unit):
                return [i, j, unit]
            j += 1
        i += 1


if __name__ == '__main__':
    cellinfo_ccle = pd.read_csv(CELLINFO_CCLE_PATH)
    tcga_drug_treatment = pd.read_csv(TCGA_DRUG_TREATMENT_PATH)
    drug_info = pd.read_csv(DRUG_INFO_PATH)
    bhk_drug_info = pd.read_csv(BHK_DRUG_INFO_PATH)

    drug_info_list = drug_info.values.tolist()
    bhk_drug_info_list = bhk_drug_info.values.tolist()

    cell_drug_list = cellinfo_ccle['treatmentid'].unique()
    tcga_drug_list = tcga_drug_treatment['pharmaceutical_therapy_drug_name'].unique()

    count_in = 0
    count_samples_in = 0
    count_not_in = 0
    count_combination = 0
    for drug in tcga_drug_list:
        found0 = search_in_list(drug, drug_info_list)
        if found0 is not None:
            print(f'found in drug_info.csv, ({found0[0]},{found0[1]}) {found0[2]}')
        found1 = search_in_list(drug, bhk_drug_info_list)
        if found1 is not None:
            print(f'found in bhk_drug_info.csv, ({found1[0]},{found1[1]}) {found1[2]}')

        if '+' in drug:
            count_combination += 1
        elif found0 is None and found1 is None:
            count_not_in += 1
        else:
            count_samples_in += tcga_drug_treatment['pharmaceutical_therapy_drug_name'].value_counts()[drug]
            count_in += 1

        # if drug in cell_drug_list:
        #     count_in += 1
        #     count_samples_in += tcga_drug_treatment['pharmaceutical_therapy_drug_name'].value_counts()[drug]
        # else:
        #     if '+' in drug:
        #         count_combination += 1
        #     else:
        #         count_not_in += 1
    print(count_in, count_not_in, count_combination, count_samples_in)
