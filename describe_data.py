import pandas as pd
from config import *

PATIENT_CELL_LINE_CSV_PATH = './final.csv'


def main():
    # df_patient_celline = pd.read_csv(PATIENT_CELL_LINE_CSV_PATH)

    df_gene_expression = pd.read_csv(CCLE_EXPRESSION_PATH)
    print(df_gene_expression.columns.values.tolist())
    # df_gene_expression_small = df_gene_expression.iloc[:, :5]
    # df_gene_expression_small.to_csv('gene_expression_small.csv', index=False)


if __name__ == '__main__':
    main()
