# Author: S.Inoue
# Date: 6/14/2021
# Updated: 6/14/2021
# Project: CTOS folfoliox project
# dataset: RioEJCA2017
# Script: Exract connected component from the network file generated through INGOR

#%%
# import module
import numpy as np
import pandas as pd
import os

#%%
def load_data(file):
    data = pd.read_table(file, sep='\t', header=0)
    print(f'[LOAD] {os.path.basename(file)}')
    print(f'# input matrix: {data.shape}')
    return data


def sample_list(data_pfsos):
    Irinotecan = ['FOLFIRI', 'FOLFIRI+BEVACIZUMAB', 'FOLFIRINOX','FOLFIRI+ERBITUX', 'FOLFIRINOX+BEVACIZUMAB', 'XELIRI+BEVACIZUMAB']
    Oxaliplatin = ['FOLFOX', 'FOLFOX+BEVACIZUMAB','FOLFIRINOX', 'FOLFIRINOX+BEVACIZUMAB']
    data_Irinotecan = data_pfsos[data_pfsos["regimen"].isin(Irinotecan)]
    data_Oxaliplatin = data_pfsos[data_pfsos["regimen"].isin(Oxaliplatin)]
    samplelist_Irinotecan = list(data_Irinotecan['Sample_name'])
    samplelist_Oxaliplatin = list(data_Oxaliplatin['Sample_name'])
    return samplelist_Irinotecan, samplelist_Oxaliplatin


def edge_selection(data_ECv, data_edges, samplelist):
    # sample selection
    l = samplelist + ['Parent', 'Child', 'name']
    data_selectedECv_ = data_ECv.loc[:,data_ECv.columns.str.contains("|".join(l))]
    # make available edgelist
    data_edges['name'] = data_edges['Parent'] + ':' + data_edges['Child']
    data_selectedECv = data_selectedECv_[data_selectedECv_['name'].isin(list(data_edges['name']))]
    print(f'# output matrix: {data_selectedECv.shape}')
    return data_selectedECv

#%%
if __name__ == '__main__':
    # load ECv data
    path_ECv = 'BayesianNetworkEstimation/CRC/ECv_Extrapolation_RioEJCA2019/ECv_extrapolation_RioEJCA2017_dataset.txt'
    data_ECv = load_data(path_ECv)

    # make sample list by each drug
    path_pfsos = 'data_RioEJCA2017/pfsos_RioEJCA2017.txt'
    data_pfsos = load_data(path_pfsos)
    samplelist_Irinotecan, samplelist_Oxaliplatin = sample_list(data_pfsos)

    #================================
    # feature extravcion by each drug
    ## Irinotecan
    path_Irinotecan = 'data/SelectedEdges_CorrCoef_rangeECv/SelectedEdges_CorrCoef0.6_rangeECv1.0_Irinotecan_10mg.txt'
    data_edges_Irinotecan = load_data(path_Irinotecan)
    # drug name
    drug_Irinotecan = path_Irinotecan.split("0_")[1].split(".")[0]
    print(f'# drug : {drug_Irinotecan}')
    # edge selection
    data_selectedECv_Irinotecan = edge_selection(data_ECv, data_edges_Irinotecan, samplelist_Irinotecan)
    # savedata
    savepath = f'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_{drug_Irinotecan}.txt'
    data_selectedECv_Irinotecan.to_csv(savepath, sep='\t', index=True)
    print(f'[SAVE] {savepath}')

    ## Oxaliplatin
    path_Oxaliplatin = 'data/SelectedEdges_CorrCoef_rangeECv/SelectedEdges_CorrCoef0.6_rangeECv1.0_Oxaliplatin_10mg.txt'
    data_edges_Oxaliplatin = load_data(path_Oxaliplatin)
    # drug name
    drug_Oxaliplatin = path_Oxaliplatin.split("0_")[1].split(".")[0]
    print(f'# drug : {drug_Oxaliplatin}')
    # edge selection
    data_selectedECv_Oxaliplatin = edge_selection(
        data_ECv, data_edges_Oxaliplatin, samplelist_Oxaliplatin)
    # savedata
    savepath = f'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_{drug_Oxaliplatin}.txt'
    data_selectedECv_Oxaliplatin.to_csv(savepath, sep='\t', index=True)
    print(f'[SAVE] {savepath}')
