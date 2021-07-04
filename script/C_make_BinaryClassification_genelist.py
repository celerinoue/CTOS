# Author: S.Inoue
# Date: 05/10/2021
# Updated: 05/10/2021
# Project: CTOS folfoliox project
# Script: To perform basic analyses

#%%
# import module
import numpy as np
import pandas as pd
import glob
import networkx as nx
import itertools
import os

#%%
def extract_Genelist_ECv(file):
    # ECv Data
    # Selected Edges (th = {CorrCoef 0.6 & ECv >= 1.0} )

    # load data
    data = pd.read_table(file, sep='\t', header=0)
    print(f'[LOAD] {os.path.basename(file)}')  # file name
    Genelist_ECv_th06 = list(set(data['Parent']) | set(data['Child']))
    print(f'# selected gene : {len(Genelist_ECv_th06)}')
    return Genelist_ECv_th06


#%%
def extract_Genelist_GeneExp(file):
    # ECv Data
    # Selected Edges (th = {CorrCoef 0.7+0.8 & ECv >= 1.0} )

    # load data
    data = pd.read_table(file, sep='\t', header=0)
    print(f'[LOAD] {os.path.basename(file)}')  # file name

    # extract gene list by every connected component
    Genelist_GeneExp = data["GeneName"].to_list()
    print(f'# selected gene : {len(Genelist_GeneExp)}')
    return Genelist_GeneExp


#%%
if __name__ == '__main__':

    # extract_Genelist_ECv
    for f in sorted(glob.glob("data/SelectedEdges_CorrCoef_rangeECv/*.txt")):
        drug = f.split("0_")[1].split(".")[0]
        print(f'# drug : {drug}')  # drug name
        if drug == 'Cetuximab_60mg':
            Genelist_ECv_th06_Cetuximab = extract_Genelist_ECv(f)
            # save matrix
            savepath = f'data_GeneList/Genelist_ECv_th06_Cetuximab.txt'
            pd.Series(Genelist_ECv_th06_Cetuximab).to_csv(
                savepath, sep='\t', header=False, index=False)
            print(f'[SAVE] {savepath}')
        elif drug == 'Irinotecan_10mg':
            Genelist_ECv_th06_Irinotecan = extract_Genelist_ECv(f)
            savepath = f'data_GeneList/Genelist_ECv_th06_Irinotecan.txt'
            pd.Series(Genelist_ECv_th06_Irinotecan).to_csv(
                savepath, sep='\t', header=False, index=False)
            print(f'[SAVE] {savepath}')
        elif drug == 'Oxaliplatin_10mg':
            Genelist_ECv_th06_Oxaliplatin = extract_Genelist_ECv(f)
            savepath = f'data_GeneList/Genelist_ECv_th06_Oxaliplatin.txt'
            pd.Series(Genelist_ECv_th06_Oxaliplatin).to_csv(
                savepath, sep='\t', header=False, index=False)
            print(f'[SAVE] {savepath}')

    # extract_Genelist_GeneExp
    for f in sorted(glob.glob("data/SelectedGeneSet_CorrCoef_rangeGeneExpression/*.txt")):
        drug = f.split("0_")[1].split(".")[0]
        print(f'# drug : {drug}')  # drug name

        if '0.7' in f:
            if drug == 'Cetuximab_60mg':
                Genelist_GeneExp_th07_Cetuximab = extract_Genelist_GeneExp(f)
                # save matrix
                savepath = f'data_GeneList/Genelist_GeneExp_th07_Cetuximab.txt'
                pd.Series(Genelist_GeneExp_th07_Cetuximab).to_csv(
                    savepath, sep='\t', header=False, index=False)
                print(f'[SAVE] {savepath}')
            elif drug == 'Irinotecan_10mg':
                Genelist_GeneExp_th07_Irinotecan = extract_Genelist_GeneExp(f)
                # save matrix
                savepath = f'data_GeneList/Genelist_GeneExp_th07_Irinotecan.txt'
                pd.Series(Genelist_GeneExp_th07_Irinotecan).to_csv(
                    savepath, sep='\t', header=False, index=False)
                print(f'[SAVE] {savepath}')
            elif drug == 'Oxaliplatin_10mg':
                Genelist_GeneExp_th07_Oxaliplatin = extract_Genelist_GeneExp(f)
                # save matrix
                savepath = f'data_GeneList/Genelist_GeneExp_th07_Oxaliplatin.txt'
                pd.Series(Genelist_GeneExp_th07_Oxaliplatin).to_csv(
                    savepath, sep='\t', header=False, index=False)
                print(f'[SAVE] {savepath}')
        if '0.8' in f:
            if drug == 'Cetuximab_60mg':
                Genelist_GeneExp_th08_Cetuximab = extract_Genelist_GeneExp(f)
                # save matrix
                savepath = f'data_GeneList/Genelist_GeneExp_th08_Cetuximab.txt'
                pd.Series(Genelist_GeneExp_th08_Cetuximab).to_csv(
                    savepath, sep='\t', header=False, index=False)
                print(f'[SAVE] {savepath}')
            elif drug == 'Irinotecan_10mg':
                Genelist_GeneExp_th08_Irinotecan = extract_Genelist_GeneExp(f)
                # save matrix
                savepath = f'data_GeneList/Genelist_GeneExp_th08_Irinotecan.txt'
                pd.Series(Genelist_GeneExp_th08_Irinotecan).to_csv(
                    savepath, sep='\t', header=False, index=False)
                print(f'[SAVE] {savepath}')
            elif drug == 'Oxaliplatin_10mg':
                Genelist_GeneExp_th08_Oxaliplatin = extract_Genelist_GeneExp(f)
                # save matrix
                savepath = f'data_GeneList/Genelist_GeneExp_th08_Oxaliplatin.txt'
                pd.Series(Genelist_GeneExp_th08_Oxaliplatin).to_csv(
                    savepath, sep='\t', header=False, index=False)
                print(f'[SAVE] {savepath}')

    # set1
    set1 = set(Genelist_ECv_th06_Cetuximab) & set(
        Genelist_GeneExp_th07_Cetuximab)
    print(f'# selected gene set1 : {len(set1)}')
    # save matrix
    savepath = f'data_GeneList/Genelist_ECv06&GeneExp07_Cetuximab.txt'
    pd.Series(list(set1)).to_csv(savepath, sep='\t', header=False, index=False)
    print(f'[SAVE] {savepath}')
    # set2
    set2 = set(Genelist_ECv_th06_Irinotecan) & set(
        Genelist_GeneExp_th07_Irinotecan)
    print(f'# selected gene set2 : {len(set2)}')
    savepath = f'data_GeneList/Genelist_ECv06&GeneExp07_Irinotecan.txt'
    pd.Series(list(set2)).to_csv(savepath, sep='\t', header=False, index=False)
    print(f'[SAVE] {savepath}')
    # set3
    set3 = set(Genelist_ECv_th06_Oxaliplatin) & set(
        Genelist_GeneExp_th07_Oxaliplatin)
    print(f'# selected gene set3 : {len(set3)}')
    savepath = f'data_GeneList/Genelist_ECv06&GeneExp07_Oxaliplatin.txt'
    pd.Series(list(set3)).to_csv(savepath, sep='\t', header=False, index=False)
    print(f'[SAVE] {savepath}')
    # set4
    set4 = set(Genelist_ECv_th06_Cetuximab) & set(
        Genelist_GeneExp_th08_Cetuximab)
    print(f'# selected gene set4 : {len(set4)}')
    savepath = f'data_GeneList/Genelist_ECv06&GeneExp08_Cetuximab.txt'
    pd.Series(list(set4)).to_csv(savepath, sep='\t', header=False, index=False)
    print(f'[SAVE] {savepath}')
    # set5
    set5 = set(Genelist_ECv_th06_Irinotecan) & set(
        Genelist_GeneExp_th08_Irinotecan)
    print(f'# selected gene set5 : {len(set5)}')
    savepath = f'data_GeneList/Genelist_ECv06&GeneExp08_Irinotecan.txt'
    pd.Series(list(set5)).to_csv(savepath, sep='\t', header=False, index=False)
    print(f'[SAVE] {savepath}')
    # set6
    set6 = set(Genelist_ECv_th06_Oxaliplatin) & set(
        Genelist_GeneExp_th08_Oxaliplatin)
    print(f'# selected gene set6: {len(set6)}')
    savepath = f'data_GeneList/Genelist_ECv06&GeneExp08_Oxaliplatin.txt'
    pd.Series(list(set6)).to_csv(savepath, sep='\t', header=False, index=False)
    print(f'[SAVE] {savepath}')
# %%
