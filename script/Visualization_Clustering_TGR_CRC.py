# Author: S.Inoue
# Date: 12/17/2020
# Updated: 03/18/2020
# Project: CTOS folfoliox project
# Script: Labeling Drug Responsiveness of 3 drugs
# Process: Idea1, Step1

# import module
import numpy as np
import pandas as pd
import pickle
import itertools
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette


def dataload(filename):
    data = pd.read_csv(filename)
    print(f'[LOAD]: {filename}')
    print(f'input matrix: {data.shape}')
    return data


def clustering(data):
    ## make list##
    drug_list = data["drug_id"].drop_duplicates().tolist()
    patient_list = data["patient_id"].drop_duplicates().tolist()
    day_list = data["day"].drop_duplicates().tolist()
    drug_index = data.drop_duplicates('drug_id').iloc[:, 1:3].reset_index(drop=True)

    l_value = []
    for l_drug, l_patient in itertools.product(drug_list, patient_list):
        l_ = data.query(f'(drug_id == {l_drug}) & (patient_id == {l_patient})').groupby(
            'day')['value'].median().reset_index(drop=True)
        l_value.append(l_)
    l_value_ = list(np.array_split(l_value, 7))
    l_gr = []
    for i in drug_list:
        gr_ = pd.DataFrame(l_value_[i-1], index=patient_list, columns=day_list)
        l_gr.append(gr_)

    ## ##
    for i in range(len(drug_index)):
        # read data
        d = l_gr[i]
        labels = list(d.index)
        # clustering (method = ward)
        z = linkage(d, method='ward')
        # plot
        plt.figure(figsize=(30, 15))
        plt.rcParams["font.family"] = 'sans-serif'
        plt.rcParams['font.size'] = 30
        plt.title(
            f'Clustering Dendrogram : {drug_index["drug_name"][i]}', fontsize=50)
        plt.xlabel('CTOS_line', fontsize=40)
        plt.ylabel('Distance', fontsize=40)
        set_link_color_palette(['red', 'blue'])  # 2 cluster
        dendrogram(z, leaf_font_size=40, color_threshold=7.,
                   labels=labels, above_threshold_color='black')
        # save & show
        savepath = f'fig/TumorGrowthRate/Clustering_TGR/fig_clustering_{drug_index["drug_name"][i]}.png'
        plt.savefig(savepath, dpi=300, format='png',
                    bbox_inches="tight")  # save
        print(f'[SAVE]: {savepath}')
    return


if __name__ == '__main__':
    # data load
    TGR_data = dataload('data_TGR/TGR.csv')
    # plot
    clustering(TGR_data)
