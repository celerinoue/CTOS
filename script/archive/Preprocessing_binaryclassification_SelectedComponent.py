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


#%%
def pickup_gene_list():
    # [LOAD] Selected ECv Edges (threshold = CorrCoef0.6) ===========
    list_gene_all = []
    for i in sorted(glob.glob("data/SelectedEdges_CorrCoef_rangeECv/*.txt")):
        # load data
        print(f'[LOAD] {os.path.basename(i)}')
        data = pd.read_table(i, sep='\t', header=0)
        # extract connected component
        G = nx.from_pandas_edgelist(data, source='Parent', target='Child', edge_attr=None, create_using=nx.MultiGraph())
        print(f'# nodes: {nx.number_of_nodes(G)}')
        print(f'# edges: {nx.number_of_edges(G)}')
        print(f'# total connected components: {nx.number_connected_components(G)}')

        list_ecv = []
        list_gene_all_ = []
        for i in sorted(nx.connected_components(G), key=len, reverse=True):
            Ge = G.subgraph(i)  # each connected graph
            if nx.number_of_nodes(Ge) >= 4:
                print(f'# component {i}')
                # list_ecv
                list_ecv.append(i)
                # list_GeneExp 発現値
                list_GeneExp = []
                for l in Ge.edges():
                    Parent_Child = str(l[0]) + "::" + str(l[1])
                    Child_Parent = str(l[1]) + "::" + str(l[0])
                    list_GeneExp.append(Parent_Child)
                    list_GeneExp.append(Child_Parent)
                list_gene_all_.append([list_ecv, list_GeneExp])
            else:
                pass
        list_gene_all.append(list_gene_all_) # len=3, (Cet=14,Irn=6,Oxa=10)
    return list_gene_all


def load_CTOS_data():
    # [LOAD] original ECv data ===========
    df_ecv_path = 'BayesianNetworkEstimation/CRC/ECv_network_result_0.15.txt'
    df_ecv_ = pd.read_table(df_ecv_path, sep='\t', header=0)
    print(f'[LOAD] {df_ecv_path}, input matrix: {df_ecv_.shape}')

    # reshape data_ecv
    df_ecv_['Parent_Child'] = [
        f"{df_ecv_['Parent'][i]}::{df_ecv_['Child'][i]}" for i in range(len(df_ecv_))]
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
    print(f'# reshaped : {df_ecv_path}')

    # [LOAD] original gene expression data ============
    df_GeneExp_path = 'BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt'
    df_GeneExp_ = pd.read_table(df_GeneExp_path, sep='\t', header=0)
    print(f'[LOAD] {df_GeneExp_path}, input matrix: {df_GeneExp_.shape}')
    df_GeneExp = df_GeneExp_.loc[:, ['GeneName',
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
    print(f'# reshaped :{df_GeneExp_path}')
    return df_ecv, df_GeneExp


def save_matrix(df_ecv, df_GeneExp, list_gene_all):
    # set the ylabels
    label_cxm60 = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1]  # 2,10
    label_irn10 = [0, 0, 0, 1, 0, 1, 0, 0, 0, 1]  # 4,6,10
    label_oxa10 = [0, 0, 1, 1, 0, 0, 1, 0, 1, 1]  # 3,4,7,9,10
    label = [label_cxm60, label_irn10, label_oxa10]  # to list
    drug = ['Cetuximab_60mg', 'Irinotecan_10mg', 'Oxaliplatin_10mg']

    # (r: drug, i: each component)
    for r in range(len(list_gene_all)): # drug
        for i in range(len(list_gene_all[r])):
            # gene_expression
            labeled_GeneExp_ = df_GeneExp[df_GeneExp['GeneName'].isin(list_gene_all[r][i][0][i])].set_index('GeneName').T
            labeled_GeneExp_["label"] = label[r]  # set ylabel
            labeled_GeneExp = labeled_GeneExp_
            #labeled_GeneExp = labeled_GeneExp_.reset_index()
            savepath = f'data_BinaryClassification/SelectedComponent/ClassificationDataSet_GeneExp_{drug[r]}_rank_{i+1}.txt'
            labeled_GeneExp.to_csv(savepath, sep='\t', index=True)
            print(f'[SAVE] {savepath}')

            # ecv : Parent_Child
            labeled_ecv_ = df_ecv[df_ecv['Parent_Child'].isin(list_gene_all[r][i][1])].set_index('Parent_Child').T
            labeled_ecv_["label"] = label[r]  # set ylabel
            #labeled_ecv = labeled_ecv_.reset_index()
            labeled_ecv = labeled_ecv_
            savepath = f'data_BinaryClassification/SelectedComponent/ClassificationDataSet_ECv_{drug[r]}_rank_{i+1}.txt'
            labeled_ecv.to_csv(savepath, sep='\t', index=True)
            print(f'[SAVE] {savepath}')
    return

if __name__ == '__main__':
    # load selected gene data
    list_gene_all = pickup_gene_list()
    # load original CTOS data
    df_ecv, df_GeneExp = load_CTOS_data()
    # reshape & save matrix
    save_matrix(df_ecv, df_GeneExp, list_gene_all)
