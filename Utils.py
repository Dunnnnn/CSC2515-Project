import copy
import random

import numpy as np
import pandas as pd
from config import *
from Record import CellineData

drug_synonyms_list = None
bhk_drug_synonyms_list = None

drug_used_by_cellines_set = None

dict_celline_ACH_XXXXXX_CCLE_synonyms = None

sampleid_ccle_synonyms_list = None
dict_ach_xxxxxx_sampleid = None

drugid_sampleid_drug_name_list = None

dict_drugid_ccle_sensitivity = None


def generate_drug_synonyms_list():
    df_drug_synonyms = pd.read_csv(DRUG_INFO_PATH)
    df_drug_synonyms['Synonyms'] = df_drug_synonyms['Synonyms'].str.split(',')
    l0 = df_drug_synonyms['Name'].values.tolist()
    l2 = df_drug_synonyms['Synonyms'].values.tolist()
    result = []
    for name0, names in zip(l0, l2):
        tmplist = [get_processed_string(name0)]
        if type(names) is list:
            for name in names:
                processed_name = get_processed_string(name)
                if processed_name not in tmplist:
                    tmplist.append(processed_name)
        result.append(tmplist)
    return result


def generate_bhk_drug_synonyms_list():
    df_bhk_drug_synonyms = pd.read_csv(BHK_DRUG_INFO_PATH)
    df_bhk_drug_synonyms['combined'] = df_bhk_drug_synonyms.apply(
        lambda x: list([x['Unnamed: 0'],
                        x['Compound..code.or.generic.name.'],
                        x['treatmentid'],
                        x['drug.name']]), axis=1)
    result = []
    for _l in df_bhk_drug_synonyms['combined'].values.tolist():
        tmp = []
        for name in _l:
            processed_name = get_processed_string(name)
            if processed_name not in tmp:
                tmp.append(processed_name)
        result.append(tmp)
    return result


def get_dictionary(df, col_name0, col_name1):
    result = df.groupby(col_name0)[col_name1].apply(list).to_dict()
    return result


def print_dic(_dic):
    for key in _dic:
        if len(_dic[key]) != 1:
            print(f'len = {len(_dic[key])}')
        print(key, _dic[key], type(_dic[key][0]))


def get_processed_string(input_string):
    tmp = []
    for c in input_string:
        if c.isalpha():
            tmp.append(c.lower())
        elif c.isdigit():
            tmp.append(c)
    return ''.join(tmp)


def find_name_in_synonyms_list(name, synonyms_list):
    for names in synonyms_list:
        for n in names:
            if name == n:
                return names
    return None


def merge_list(list0, list1):
    result = copy.copy(list0)
    result.extend(x for x in list1 if x not in result)
    return result


def get_drug_synonyms_list(drug_name):
    global drug_synonyms_list, bhk_drug_synonyms_list
    if drug_synonyms_list is None or bhk_drug_synonyms_list is None:
        drug_synonyms_list = generate_drug_synonyms_list()
        bhk_drug_synonyms_list = generate_bhk_drug_synonyms_list()

    result0 = find_name_in_synonyms_list(drug_name, drug_synonyms_list)
    if result0 is None:
        result0 = [drug_name]

    for name in result0:
        result1 = find_name_in_synonyms_list(name, bhk_drug_synonyms_list)
        if result1 is not None:
            result0 = merge_list(result0, result1)

    return sorted(result0)


def str_list_2_key(str_list):
    return ','.join(str_list)


def key_2_str_list(key):
    return key.split(',')


def is_in_drug_used_by_cellines_set(drug_name):
    global drug_used_by_cellines_set
    if drug_used_by_cellines_set is None:
        drug_used_by_cellines_set = set([])
        df_cellinfo_ccle = pd.read_csv(CELL_INFO_CCLE_PATH)
        for treatment in df_cellinfo_ccle['treatmentid'].unique():
            drug_used_by_cellines_set.add(get_processed_string(treatment))
    return drug_name in drug_used_by_cellines_set


def get_CCLE_synonyms_with_ach_xxxxxx(ach_xxxxxx):
    global dict_celline_ACH_XXXXXX_CCLE_synonyms
    if dict_celline_ACH_XXXXXX_CCLE_synonyms is None:
        dict_celline_ACH_XXXXXX_CCLE_synonyms = {}
        df_sample_info = pd.read_csv(SAMPLE_INFO_PATH)
        df_sample_info['combined'] = df_sample_info.apply(
            lambda x: list([x['stripped_cell_line_name'],
                            x['CCLE_Name']]), axis=1)
        ach_xxxxxx_list = df_sample_info['DepMap_ID'].values.tolist()
        ccle_synonyms_list = df_sample_info['combined'].values.tolist()
        for ach_id, ccle_synonyms in zip(ach_xxxxxx_list, ccle_synonyms_list):
            tmp_list = []
            for synonyms in ccle_synonyms:
                if type(synonyms) is str:
                    tmp_list.append(get_processed_string(synonyms))
            dict_celline_ACH_XXXXXX_CCLE_synonyms[ach_id] = tmp_list

    return dict_celline_ACH_XXXXXX_CCLE_synonyms.get(ach_xxxxxx, None)


def has_match(list0, list1):
    for s0 in list0:
        for s1 in list1:
            if s0 == s1:
                return True
    return False


def get_sampleid_with_ach_xxxxxx(ach_xxxxxx):
    global sampleid_ccle_synonyms_list, dict_ach_xxxxxx_sampleid

    # cached
    if dict_ach_xxxxxx_sampleid is not None:
        matched_sampleid = dict_ach_xxxxxx_sampleid.get(ach_xxxxxx)
        if matched_sampleid is not None:
            return matched_sampleid

    ccle_synonyms0 = get_CCLE_synonyms_with_ach_xxxxxx(ach_xxxxxx)
    if ccle_synonyms0 is None:
        return None

    if sampleid_ccle_synonyms_list is None:
        sampleid_ccle_synonyms_list = []
        df_bhk_cell_info = pd.read_csv(BHK_CELL_INFO_PATH)
        sampleid_list = df_bhk_cell_info['sampleid']
        ccle_name_list = df_bhk_cell_info['CCLE.name']
        ccle_pri_name_list = df_bhk_cell_info['Cell.line.primary.name']
        for _sampleid, _ccle_name, _ccle_pri_name in zip(sampleid_list, ccle_name_list, ccle_pri_name_list):
            if type(_sampleid) is str:
                tmp = []
                if type(_ccle_name) is str:
                    tmp.append(get_processed_string(_ccle_name))
                if type(_ccle_pri_name) is str:
                    tmp.append(get_processed_string(_ccle_pri_name))
                if len(tmp) > 0:
                    sampleid_ccle_synonyms_list.append((get_processed_string(_sampleid), tmp))

    matched_sampleid = None
    for sampleid, ccle_synonyms_list in sampleid_ccle_synonyms_list:
        if has_match(ccle_synonyms0, ccle_synonyms_list):
            matched_sampleid = sampleid
            break

    if matched_sampleid is None:
        return None

    if dict_ach_xxxxxx_sampleid is None:
        dict_ach_xxxxxx_sampleid = {}
    dict_ach_xxxxxx_sampleid[ach_xxxxxx] = matched_sampleid
    return matched_sampleid


def find_drugid_with_sampleid_n_drug_name(sampleid_input, drug_name):
    global drugid_sampleid_drug_name_list
    if drugid_sampleid_drug_name_list is None:
        drugid_sampleid_drug_name_list = []
        df_cellinfo_ccle = pd.read_csv(CELL_INFO_CCLE_PATH)
        drugid_list = df_cellinfo_ccle['Unnamed: 0'].values.tolist()
        sampleid_list = df_cellinfo_ccle['sampleid'].values.tolist()
        treatmentid_list = df_cellinfo_ccle['treatmentid'].values.tolist()
        for drugid, sampleid, treatmentid in zip(drugid_list, sampleid_list, treatmentid_list):
            drugid_sampleid_drug_name_list.append([drugid,
                                                   get_processed_string(sampleid), get_processed_string(treatmentid)])

    for drugid, sampleid, treatmentid in drugid_sampleid_drug_name_list:
        if sampleid == sampleid_input and treatmentid == drug_name:
            return drugid

    return None


def get_celline_sensitivity_with_drugid(drugid):
    global dict_drugid_ccle_sensitivity
    if dict_drugid_ccle_sensitivity is None:
        df_aucs_ccle = pd.read_csv(AUCS_CCLE_PATH)
        dict_drugid_ccle_sensitivity = df_aucs_ccle.set_index('Unnamed: 0').to_dict('index')

    return dict_drugid_ccle_sensitivity.get(drugid, None)


def find_cellines_with_drug(key_drug_synonyms, cell_line_ach_xxxxxx_list):
    result = []
    # cell_line_ACH_XXXXXX_list = pd.read_csv(TUMOR_CL_CORRELATION_PATH, index_col=0, nrows=0).columns.tolist()
    # key_drug_synonyms = 'abraxane,bms18133901,onxol,paclitaxel,paxene,praxel,taxol'

    for drug in key_2_str_list(key_drug_synonyms):
        if not is_in_drug_used_by_cellines_set(drug):
            continue
        for index, ach_xxxxxx in enumerate(cell_line_ach_xxxxxx_list):
            sampleid = get_sampleid_with_ach_xxxxxx(ach_xxxxxx)
            if sampleid is None:
                continue

            matched_drugid = find_drugid_with_sampleid_n_drug_name(sampleid, drug)
            if matched_drugid is None:
                continue

            celline_sensitivity = get_celline_sensitivity_with_drugid(matched_drugid)
            if celline_sensitivity is None:
                continue

            contains_nan = False
            for k in celline_sensitivity:
                if pd.isna(celline_sensitivity[k]):
                    contains_nan = True

            if contains_nan:
                continue

            result.append(CellineData(index=index, drug_name=drug, sampleid=sampleid, sensitivity=celline_sensitivity))

    return result


def show_cor_distribution(record_list):
    pass


def get_merged_class_name(class_str):
    if class_str == 'Stable Disease' or class_str == 'Persistent Disease':
        return 'Stable Disease/Persistent Disease'

    if class_str == 'Partial Response' or class_str == 'Partial Remission/Response':
        return 'Partial Remission/Response'

    if class_str == 'Complete Remission/Response' or class_str == 'Complete Response':
        return 'Complete Remission/Response'

    return class_str


def down_sampling(record_list):
    dict_samples = {}
    for record in record_list:
        class_name = record.recist_score
        class_record_list = dict_samples.get(class_name, None)
        if class_record_list is None:
            dict_samples[class_name] = []
        dict_samples[class_name].append(record)

    min_samples = 1000000
    for k in dict_samples:
        min_samples = min(len(dict_samples[k]), min_samples)

    random.seed(2021)
    result = []
    for k in dict_samples:
        if len(dict_samples[k]) == min_samples:
            result.extend(dict_samples[k])
        else:
            result.extend(random.sample(dict_samples[k], min_samples))
    return result
