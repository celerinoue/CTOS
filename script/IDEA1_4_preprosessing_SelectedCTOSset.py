# Author: S.Inoue
# Date: 01/29/2021
# Updated: 01/29/2021
# Project: CTOS folfoliox project


# import module
import numpy as np
import pandas as pd
import glob
import itertools
import os

# pickup gene list from selected edge file
def pickup_gene_list():

    # ECv data
    genelist_ecv_pair_th06 = []
    genelist_ecv_th06 = []
    for i in sorted(glob.glob("data/SelectedEdges_CorrCoef_rangeECv/*.txt")):
        print(f'[LOAD] {os.path.basename(i)}')
        df = pd.read_table(i, sep='\t', header=0)
        # list of gene_name of ECv pair
        df['Parent_Child'] = df['Parent'] + '::' + df['Child']  # pair list
        genelist_ecv_pair_th06.append(list(df['Parent_Child']))
        # list of gene_name of ECv union
        list_ecv = list(set(df['Parent']) | set(df['Child']))
        genelist_ecv_th06.append(list(list_ecv))

    # gene expression data
    genelist_exp_th07 = []
    genelist_exp_th08 = []
    for i in sorted(glob.glob("data/SelectedGeneSet_CorrCoef_rangeGeneExpression/*.txt")):
        print(f'[LOAD] {os.path.basename(i)}')
        df = pd.read_table(i, sep='\t', header=0)
        if '0.7' in i:
            # list of gene_name of th0.7
            genelist_exp_th07.append(list(df["GeneName"]))
        elif '0.8' in i:
            # list of gene_name of th0.8
            genelist_exp_th08.append(list(df["GeneName"]))

    gene_list_all = pd.DataFrame([genelist_ecv_pair_th06,
                                    genelist_ecv_th06,
                                    genelist_exp_th07,
                                    genelist_exp_th08],
                                index=["CorrCoef0.6_rangeECv1.0",
                                        "CorrCoef0.6_rangeECv1.0_forGeneExp",
                                        "CorrCoef0.7_rangeGeneExp1.0",
                                        "CorrCoef0.8_rangeGeneExp1.0"],
                                columns=["Cetuximab_60mg",
                                        "Irinotecan_10mg",
                                        "Oxaliplatin_10mg"])

    return gene_list_all


def load_CTOS_data():
    # [LOAD] original ECv data ===========
    df_ecv_path = 'data/ECv_network_result_0.15.txt'
    df_ecv_ = pd.read_table(df_ecv_path, sep='\t', header=0)
    print(f'[LOAD] {df_ecv_path}, input matrix: {df_ecv_.shape}')
    # reshape data_ecv
    df_ecv_['Parent_Child'] = df_ecv_['Parent'] + '::' + df_ecv_['Child']
    df_ecv = df_ecv_.loc[:, ['Parent_Child',
                             'ECv:C97-float:8',  # CTOS_line = 1
                             'ECv:C166-float:21',  # CTOS_line = 2
                             'ECv:C86-float:17',  # CTOS_line = 3
                             'ECv:C111-foat:18',  # CTOS_line = 4
                             'ECv:C45-float:5',  # CTOS_line = 5
                             'ECv:C48-float:6',  # CTOS_line = 6
                             'ECv:C138-float:20',  # CTOS_line = 7
                             'ECv:CB3-float:22',  # CTOS_line = 8
                             'ECv:C75-float:7',  # CTOS_line = 9
                             'ECv:C132-float:19'  # CTOS_line = 10
                             ]]
    print(f'[INFO] reshaped :{df_ecv_path}')

    # [LOAD] original gene expression data ============
    df_exp_path = 'data/CRC_dataset.txt'
    df_exp_ = pd.read_table(df_exp_path, sep='\t', header=0)
    print(f'[LOAD] {df_exp_path}, input matrix: {df_exp_.shape}')
    df_exp = df_exp_.loc[:, ['GeneName',
                             'C97-float',  # CTOS_line = 1
                             'C166-float',  # CTOS_line = 2
                             'C86-float',  # CTOS_line = 3
                             'C111-foat',  # CTOS_line = 4
                             'C45-float',  # CTOS_line = 5
                             'C48-float',  # CTOS_line = 6
                             'C138-float',  # CTOS_line = 7
                             'CB3-float',  # CTOS_line = 8
                             'C75-float',  # CTOS_line = 9
                             'C132-float'  # CTOS_line = 10
                             ]]
    print(f'[INFO] reshaped :{df_exp_path}')

    return df_ecv, df_exp


def save_matrix(df_ecv, df_exp, gene_list_all):
    # set the ylabels
    label_cxm60 = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1]  # 2,10
    label_irn60 = [0, 0, 0, 1, 0, 1, 0, 0, 0, 1]  # 4,6,10
    label_oxa60 = [0, 0, 1, 1, 0, 0, 1, 0, 1, 1]  # 3,4,7,9,10
    label = [label_cxm60, label_irn60, label_oxa60]  # to list

    # (r: CTOS name, i: drug name)
    for r, i in itertools.product(range(len(gene_list_all.index)), range(len(gene_list_all.columns))):
        # set the ylabels
        label_cxm60 = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1]  # 2,10
        label_irn10 = [0, 0, 0, 1, 0, 1, 0, 0, 0, 1]  # 4,6,10
        label_oxa10 = [0, 0, 1, 1, 0, 0, 1, 0, 1, 1]  # 3,4,7,9,10
        label = [label_cxm60, label_irn10, label_oxa10]  # to list
        # extract CTOS data
        if r == 0:
            selected_matrix_ = df_ecv[df_ecv['Parent_Child'].isin(gene_list_all.iloc[0, i])].set_index('Parent_Child').T
            selected_matrix_["label"] = label[i]  # set ylabel
            selected_matrix = selected_matrix_.reset_index()
        else:
            selected_matrix_ = df_exp[df_exp['GeneName'].isin(gene_list_all.iloc[r, i])].set_index('GeneName').T
            selected_matrix_["label"] = label[i]  # set ylabel
            selected_matrix = selected_matrix_.reset_index()
        # save matrix
        savepath = f'result/txt/IDEA1_4/SelectedCTOSset/SelectedCTOSset_{gene_list_all.index[r]}_{gene_list_all.columns[i]}.txt'
        selected_matrix.to_csv(savepath, sep='\t', index=False)
        print(f'[SAVE] {savepath}')
    return

if __name__ == '__main__':
    # load selected gene data
    gene_list_all = pickup_gene_list()
    # load original CTOS data
    df_ecv, df_exp = load_CTOS_data()
    # reshape & save matrix
    save_matrix(df_ecv, df_exp, gene_list_all)
