# Author: yoshi and S.Inoue
# Date: 2/13/2021
# Updated: 4/18/2021
# Project: CTOS folfoli folfox
# Script: Exract connected component from the network file generated through INGOR

#%%
# import module
import numpy as np
import pandas as pd
import glob
import networkx as nx
import itertools
import os


def load_GeneExp_data(file):
    data_ = pd.read_table(file, sep='\t', header=0)
    print(f'[LOAD] {os.path.basename(file)}')
    print(f'# input matrix: {data_.shape}')
    # reshape data
    data = data_.loc[:, ['GeneName',
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
    return data


def extract_Genelist(file):
    # ECv Data
    # Selected Edges (th = {CorrCoef 0.7+0.8 & ECv >= 1.0} )

    # load data
    data = pd.read_table(file, sep='\t', header=0)
    print(f'[LOAD] {os.path.basename(file)}')  # file name

    # extract gene list by every connected component
    Genelist_GeneExp = data["GeneName"].to_list()
    print(f'# selected gene : {len(Genelist_GeneExp)}')

    return Genelist_GeneExp


def set_ylabel():
    # drug list
    drug_list = ['Cetuximab_60mg', 'Irinotecan_10mg', 'Oxaliplatin_10mg']
    # set ylabels for classification
    ylabel_Cetuximab_60mg = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1]  # 2,10
    ylabel_Irinotecan_10mg = [0, 0, 0, 1, 0, 1, 0, 0, 0, 1]  # 4,6,10
    ylabel_Oxaliplatin_10mg = [0, 0, 1, 1, 0, 0, 1, 0, 1, 1]  # 3,4,7,9,10

    ylabel = pd.DataFrame(
        [[ylabel_Cetuximab_60mg], [ylabel_Irinotecan_10mg], [ylabel_Oxaliplatin_10mg]],
        index=drug_list)  # to list
    return drug_list, ylabel


def makeinputdataset_ECv(genelist, data, label):
    labeled_ECv = data[data['Parent_Child'].isin(
        genelist)].set_index('Parent_Child').T
    labeled_ECv["label"] = ylabel  # set ylabel
    return labeled_ECv


def makeinputdataset_GeneExp(genelist, data, ylabel):
    # create a matrix for each component ( column = gene_list, row = LINE 1~10)
    # gene_expression ( column = gene_list, row = LINE 1~10)
    labeled_GeneExp = data[data['GeneName'].isin(
        genelist)].set_index('GeneName').T
    labeled_GeneExp["label"] = ylabel  # set ylabel

    return labeled_GeneExp


#%%
if __name__ == '__main__':

    # load GeneExp data
    path_GeneExp = 'BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt'
    data_GeneExp = load_GeneExp_data(path_GeneExp)

    # reshape & save matrix ===================

    # set ylabels for classification
    drug_list, ylabel_index = set_ylabel()

    # make dataset for binary classification
    print('[INFO] extract list of gene names =======================================')
    for f in sorted(glob.glob("data/SelectedGeneSet_CorrCoef_rangeGeneExpression/*.txt")):
        print('#=======================================')
        # set ylabel
        drug = f.split("0_")[1].split(".")[0]
        print(f'# drug : {drug}')  # drug name
        ylabel = ylabel_index.loc[drug][0]
        # extract list of gene names
        Genelist_GeneExp = extract_Genelist(f)

        # create matrix for each component ======================

        # GeneExp th=0.7
        if '0.7' in f:
            labeled_GeneExp = makeinputdataset_GeneExp(Genelist_GeneExp, data_GeneExp, ylabel)
            # save matrix
            savepath = f'data/data_BinaryClassification/ClassificationDataSet_GeneExp_th07/ClassificationDataSet_GeneExp_th07_{drug}_allgene.txt'
            labeled_GeneExp.to_csv(savepath, sep='\t', index=True)
            print(f'[SAVE] {savepath}')

        # GeneExp th=0.8
        elif '0.8' in f:
            labeled_ECv = makeinputdataset_GeneExp(Genelist_GeneExp, data_GeneExp, ylabel)
            # save matrix
            savepath = f'data/data_BinaryClassification/ClassificationDataSet_GeneExp_th08/ClassificationDataSet_GeneExp_th08_{drug}_allgene.txt'
            labeled_ECv.to_csv(savepath, sep='\t', index=True)
            print(f'[SAVE] {savepath}')
