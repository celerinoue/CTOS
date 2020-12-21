# Author: S.Inoue
# Date: 12/17/2020
# Updated: 12/22/2020
# Project: CTOS folfoliox project
# Script: To perform basic analyses

# import module
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette

def data_load(file, file2):
    with open(file, mode="rb") as f:
        data = pickle.load(f)
    print(f'[LOAD]: {file}')

    data2 = pd.read_csv(file2, index_col=0)
    print(f'[LOAD]: {file2}')
    #print(f'input matrix: {data.shape}')
    return data, data2


def clustering(data, drug_index):
    for i in range(len(drug_index)):
        # read data
        d = l_gr[i]
        labels = list(d.index)
        # clustering (method = ward)
        z = linkage(d, method='ward')
        # plot
        plt.figure(figsize=(30, 15))
        plt.rcParams["font.family"] = 'Times New Roman'
        plt.rcParams['font.size'] = 30
        plt.title(
            f'Clustering Dendrogram : {drug_index["drug_name"][i]}', fontsize=50)
        plt.xlabel('CTOS_line', fontsize=40)
        plt.ylabel('Distance', fontsize=40)
        set_link_color_palette(['red', 'blue'])  # 2 cluster
        dendrogram(z, leaf_font_size=40, color_threshold=7.,
                   labels=labels, above_threshold_color='black')
        # save & show
        savepath = f'./result/fig/IDEA1_1/clustering_TumorGrowthRate/fig_clustering_{drug_index["drug_name"][i]}.png'
        plt.savefig(savepath, dpi=300, format='png',
                    bbox_inches="tight")  # save
        print(f'[SAVE]: {savepath}')
    return



if __name__ == '__main__':
    # data load
    l_gr, drug_index = data_load(
        './result/txt/IDEA1_1/tumor_growth_rate.pickle',
        './result/txt/IDEA1_1/drug_index.csv')

    # plot
    clustering(l_gr, drug_index)
