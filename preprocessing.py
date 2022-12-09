import copy
import numpy as np
import pandas as pd

from Record import PatientRecord, CellineData, FinalRecord
from config import *
from Utils import *

if __name__ == '__main__':
    df_tumor_cl_correlation = pd.read_csv(TUMOR_CL_CORRELATION_PATH)
    print(df_tumor_cl_correlation.shape)
    # select rows with TCGA in it
    df_tumor_cl_correlation = df_tumor_cl_correlation.loc[df_tumor_cl_correlation['Unnamed: 0'].str.contains('TCGA')]
    df_tumor_cl_correlation = df_tumor_cl_correlation.sort_values(by='Unnamed: 0')
    print(df_tumor_cl_correlation.shape)

    # Convert TCGA-XX-XXXX-0X to TCGA-XX-XXXX
    patient_record_list = []
    tmp = []
    for row in df_tumor_cl_correlation.values.tolist():
        tcga_id = row[0][:-3]
        if tcga_id not in tmp:
            tmp.append(tcga_id)
            patient_record_list.append(PatientRecord(patient_id=tcga_id, correlation=row[1:]))
    print(f'{len(patient_record_list)} different TCGA patients')

    # find drug used and RECIST score
    df_patient_recist = pd.read_csv(DRUG_RESPONSE_PATH)
    dict_patient_recist = get_dictionary(df_patient_recist, 'bcr_patient_barcode', 'treatment_outcome_at_tcga_followup')
    df_patient_drug = pd.read_csv(DRUG_TREATMENT_PATH)
    dict_patient_drug = get_dictionary(df_patient_drug, 'bcr_patient_barcode', 'pharmaceutical_therapy_drug_name')
    for patient_record in patient_record_list:
        drug_recist_score = dict_patient_recist.get(patient_record.patient_id, [None])[0]
        drug_treatment = dict_patient_drug.get(patient_record.patient_id, [None])[0]
        # filter out patient not found or with multiple drugs
        if drug_recist_score is not None and drug_treatment is not None and '+' not in drug_treatment:
            processed_drug_name = get_processed_string(drug_treatment)
            patient_record.drug_synonyms_list = get_drug_synonyms_list(processed_drug_name)
            patient_record.recist_score = drug_recist_score
    patient_record_list = [r for r in patient_record_list if r.recist_score is not None]
    print(f'found {len(patient_record_list)} patients with TCGA treatment record')

    # find cell lines with drugs used by patients
    # convert drug list to string to remove duplicates
    key_drug_synonyms_list = []
    for patient_record in patient_record_list:
        key_drug_synonyms = str_list_2_key(patient_record.drug_synonyms_list)
        if key_drug_synonyms not in key_drug_synonyms_list:
            key_drug_synonyms_list.append(key_drug_synonyms)
    print(f'{len(key_drug_synonyms_list)} different drugs')

    # for every drug, find cell lines, build a dict
    dict_drug_celline_sensitivity = {}
    cell_line_ach_xxxxxx_list = df_tumor_cl_correlation.columns.values[1:]
    for key_drug_synonyms in key_drug_synonyms_list:
        celline_sensitivity_list = find_cellines_with_drug(key_drug_synonyms, cell_line_ach_xxxxxx_list)
        if len(celline_sensitivity_list) > 0:
            dict_drug_celline_sensitivity[key_drug_synonyms] = celline_sensitivity_list
    print(f'{len(dict_drug_celline_sensitivity)} drugs are found in cell lines')

    # match patient and cell lines with the same drug
    final_records = []
    patient_count = 0
    for patient_record in patient_record_list:
        key_drug_synonyms = str_list_2_key(patient_record.drug_synonyms_list)
        celline_sensitivity_list = dict_drug_celline_sensitivity.get(key_drug_synonyms, None)
        if celline_sensitivity_list is not None:
            patient_count += 1
            for celline_data in celline_sensitivity_list:
                celline_name = cell_line_ach_xxxxxx_list[celline_data.index]
                cor = patient_record.correlation[celline_data.index]
                merged_class_name = get_merged_class_name(patient_record.recist_score)
                final_record = FinalRecord(patient_id=patient_record.patient_id,
                                           drug_name=celline_data.drug_name,
                                           celline=celline_name,
                                           correlation=cor,
                                           sensitivity=celline_data.sensitivity,
                                           recist_score=merged_class_name)
                final_records.append(final_record)
    print(f'{patient_count} patients, {len(final_records)} records')

    final_records = down_sampling(final_records)

    file_name = 'final.csv'
    df_final = pd.DataFrame([vars(r) for r in final_records])
    df_final.to_csv(file_name, index=False)
    print(f'save to {file_name}')

    df_gene_expression = pd.read_csv(CCLE_EXPRESSION_PATH)
    print(df_gene_expression.shape)
    ach_xxxxxx_list = df_final['celline'].unique().tolist()
    df_gene_expression = df_gene_expression.loc[df_gene_expression['Unnamed: 0'].isin(ach_xxxxxx_list)]
    print(len(ach_xxxxxx_list), df_gene_expression.shape)
    df_gene_expression.to_csv('gene_expression_we_need.csv', index=False)

    # df_final_merged = df_final.merge(df_gene_expression, left_on='celline', right_on='Unnamed: 0')
    # df_final_merged = df_final_merged.drop('Unnamed: 0', axis=1)
    # file_name = 'final_merged.csv'
    # df_final_merged.to_csv(file_name, index=False)
    # print(f'save merged to {file_name}')
    #
    # print(df_final_merged.shape)
    # print(df_final_merged['recist_score'].value_counts())
