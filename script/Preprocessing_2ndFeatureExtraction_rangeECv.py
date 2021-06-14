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


def load_ECv_data(file):
    data_ = pd.read_table(file, sep='\t', header=0)
    print(f'[LOAD] {os.path.basename(file)}')
    print(f'# input matrix: {data_.shape}')
    # reshape data
    data_['Parent_Child'] = [
        f"{data_['Parent'][i]}::{data_['Child'][i]}" for i in range(len(data_))]
    data = data_.loc[:, ['Parent_Child',
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
    return data


def extract_Genelist(file):
    # ECv Data
    # Selected Edges (th = {CorrCoef 0.6 & ECv >= 1.0} )

    # load data
    data = pd.read_table(file, sep='\t', header=0)
    print(f'[LOAD] {os.path.basename(file)}') # file name

    # count number of noeds, edges, and connected components
    G = nx.from_pandas_edgelist(
        data, source='Parent', target='Child', edge_attr=None, create_using=nx.MultiGraph())
    print(f'# nodes: {nx.number_of_nodes(G)}')
    print(f'# edges: {nx.number_of_edges(G)}')
    print(f'# total connected components: {nx.number_connected_components(G)}')


    # extract gene list by every connected component

    # Genelist_GeneExp_th06_all
    Genelist_GeneExp_th06_all = list(set(data['Parent']) | set(data['Child']))
    print(f'# gene (GeneExp, th06, all component) : {len(Genelist_GeneExp_th06_all)}')

    # Genelist_ECv_th06_all
    data['Parent_Child'] = data['Parent'] + '::' + data['Child']  # pair list
    data['Child_Parent'] = data['Child'] + '::' + data['Parent']
    Genelist_ECv_th06_all = list(set(data['Parent_Child']) | set(data['Child_Parent']))
    print(f'# gene (ECv, th06, all component) : {len(Genelist_ECv_th06_all)}')

    # Genelist_GeneExp_th06, Genelist_ECv_th06
    Genelist_GeneExp_th06 = []
    Genelist_ECv_th06 = []
    for i in sorted(nx.connected_components(G), key=len, reverse=True):
        Ge = G.subgraph(i)  # each connected component
        if nx.number_of_nodes(Ge) >= 4:  # Limit to 4 or bigger components
            print(f'    # component {i}')
            # list_GeneExp
            Genelist_GeneExp_th06.append(i)
            # list_EXv
            list_ = []
            for l in Ge.edges():  # directed graph
                Parent_Child = str(l[0]) + "::" + str(l[1])
                Child_Parent = str(l[1]) + "::" + str(l[0])
                list_.append(Parent_Child)
                list_.append(Child_Parent)
            Genelist_ECv_th06.append(list_)
        else:
            pass
    print(f'# selected connected components: {len(Genelist_GeneExp_th06)}')
    #print(f'# Gene list of ECv    : {len(Genelist_ECv_th06)}')
    #print(f'# Gene list of GeneExp: {len(Genelist_GeneExp_th06)}')

    return Genelist_GeneExp_th06, Genelist_ECv_th06, Genelist_GeneExp_th06_all, Genelist_ECv_th06_all


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
    labeled_ECv = data[data['Parent_Child'].isin(genelist)].set_index('Parent_Child').T
    labeled_ECv["label"] = ylabel  # set ylabel
    return labeled_ECv


def makeinputdataset_GeneExp(genelist, data, ylabel):
    # create a matrix for each component ( column = gene_list, row = LINE 1~10)
    # gene_expression ( column = gene_list, row = LINE 1~10)
    labeled_GeneExp = data[data['GeneName'].isin(genelist)].set_index('GeneName').T
    labeled_GeneExp["label"] = ylabel  # set ylabel
    return labeled_GeneExp


if __name__ == '__main__':

    # load GeneExp data
    path_GeneExp = 'BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt'
    data_GeneExp = load_GeneExp_data(path_GeneExp)

    # load ECv data
    path_ECv = 'BayesianNetworkEstimation/CRC/ECv_network_result_0.15.txt'
    data_ECv = load_ECv_data(path_ECv)


    # reshape & save matrix ===================

    # set ylabels for classification
    drug_list, ylabel_index = set_ylabel()

    # make dataset for binary classification
    print('[INFO] extract list of gene names =======================================')
    for f in sorted(glob.glob("data/SelectedEdges_CorrCoef_rangeECv/*.txt")):
        print('#=======================================')

        # set ylabel
        drug = f.split("0_")[1].split(".")[0]
        print(f'# drug : {drug}')  # drug name
        ylabel = ylabel_index.loc[drug][0]
        # extract list of gene names
        Genelist_GeneExp_th06, Genelist_ECv_th06, Genelist_GeneExp_th06_all, Genelist_ECv_th06_all = extract_Genelist(f)

        # create matrix for each component ======================

        # GeneExp th=0.6
        for i in range(len(Genelist_GeneExp_th06)):
            labeled_GeneExp = makeinputdataset_GeneExp(Genelist_GeneExp_th06[i], data_GeneExp, ylabel)
            # save matrix
            savepath = f'data_2ndFeatureExtractedDataset/ClassificationDataSet_GeneExp_rangeECv_th06/ClassificationDataSet_GeneExp_rangeECv_th06_{drug}_rank_{i+1}.txt'
            labeled_GeneExp.to_csv(savepath, sep='\t', index=True)
            print(f'[SAVE] {savepath}')

        # ECv th=0.6
        for i in range(len(Genelist_ECv_th06)):
            labeled_ECv = makeinputdataset_ECv(Genelist_ECv_th06[i], data_ECv, ylabel)
            # save matrix
            savepath = f'data_2ndFeatureExtractedDataset/ClassificationDataSet_ECv_th06/ClassificationDataSet_ECv_th06_{drug}_rank_{i+1}.txt'
            labeled_ECv.to_csv(savepath, sep='\t', index=True)
            print(f'[SAVE] {savepath}')

        # GeneExp th=0.6 all
        labeled_GeneExp_all = makeinputdataset_GeneExp(Genelist_GeneExp_th06_all, data_GeneExp, ylabel)
        # save matrix
        savepath = f'data_2ndFeatureExtractedDataset/ClassificationDataSet_GeneExp_rangeECv_th06/ClassificationDataSet_GeneExp_rangeECv_th06_{drug}_allgene.txt'
        labeled_GeneExp_all.to_csv(savepath, sep='\t', index=True)
        print(f'[SAVE] {savepath}')

        # ECv th=0.6 all
        labeled_ECv_all = makeinputdataset_ECv(Genelist_ECv_th06_all, data_ECv, ylabel)
        # save matrix
        savepath = f'data_2ndFeatureExtractedDataset/ClassificationDataSet_ECv_th06/ClassificationDataSet_ECv_th06_{drug}_allgene.txt'
        labeled_ECv_all.to_csv(savepath, sep='\t', index=True)
        print(f'[SAVE] {savepath}')

# %%
