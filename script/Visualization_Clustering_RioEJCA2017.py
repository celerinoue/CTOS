# Author: S.Inoue
# Date: 6/14/2021
# Updated: 6/14/2021
# Project: CTOS folfoliox project
# dataset: RioEJCA2017
# Script: clustering for SurvivalAnalysis

#%%
# import module
from scipy.cluster.hierarchy import linkage, dendrogram
import numpy as np
import pandas as pd
import os
import itertools
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette

#%%
def load_data(file):
    data = pd.read_table(file, sep='\t', header=0, index_col=0)
    print(f'[LOAD] {os.path.basename(file)}')
    print(f'# input matrix: {data.shape}')
    return data
#%%
def reshape_inputmatrix(data_):
    data = data_.set_index('name').drop(columns=['Parent', 'Child']).T # (sample * edge)
    return data

#%%
def clustering(data):
    labels = list(input_data.index)
    result = linkage(data.iloc[:, :],
                    #metric = 'braycurtis',
                    #metric = 'canberra',
                    #metric = 'chebyshev',
                    #metric = 'cityblock',
                    metric='correlation',
                    #metric = 'cosine',
                    #metric = 'euclidean',
                    #metric = 'hamming',
                    #metric = 'jaccard',
                    #method= 'single')
                    method='average')
                    #method= 'complete')
                    #method='weighted')
    return result, labels


#%%
def plot(result,labels, drug, savepath):
    plt.figure(figsize=(10, 10))
    dendrogram(result, orientation='right', labels=labels, color_threshold=0.01)
    plt.title(f"Dendrogram of RioEJCA2017 SelectedSamples [drug = {drug}]")
    plt.xlabel("Threshold")
    #plt.grid()
    #plt.show()
    plt.savefig(savepath, dpi=100, format='png', bbox_inches="tight")  # save
    print(f'[SAVE]: {savepath}')
    return


#%%
def draw_threshold_dependency(result, savepath):
    n_clusters = len(result)
    n_samples = len(result)
    df1 = pd.DataFrame(result)
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    for i in range(len(result) - 1):
        n1 = int(result[i][0])
        n2 = int(result[i][1])
        val = result[i][2]
        n_clusters -= 1
        x1.append(val)
        x2.append(val)
        y1.append(n_clusters)
        y2.append(float(n_samples) / float(n_clusters))

    plt.subplot(2, 1, 1)
    plt.plot(x1, y1, 'yo-')
    plt.title('Threshold dependency of hierarchical clustering')
    plt.ylabel('Num of clusters')
    plt.subplot(2, 1, 2)
    plt.plot(x2, y2, 'ro-')
    plt.xlabel('Threshold')
    plt.ylabel('Ave cluster size')
    #plt.show()
    plt.savefig(savepath, dpi=100, format='png', bbox_inches="tight")  # save
    print(f'[SAVE]: {savepath}')
    return


def get_cluster_by_number(result, number):
    output_clusters = []
    x_result, y_result = result.shape
    n_clusters = x_result + 1
    cluster_id = x_result + 1
    father_of = {}
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    for i in range(len(result) - 1):
        n1 = int(result[i][0])
        n2 = int(result[i][1])
        val = result[i][2]
        n_clusters -= 1
        if n_clusters >= number:
            father_of[n1] = cluster_id
            father_of[n2] = cluster_id

        cluster_id += 1

    cluster_dict = {}
    for n in range(x_result + 1):
        if n not in father_of:
            output_clusters.append([n])
            continue

        n2 = n
        m = False
        while n2 in father_of:
            m = father_of[n2]
            #print [n2, m]
            n2 = m

        if m not in cluster_dict:
            cluster_dict.update({m: []})
        cluster_dict[m].append(n)

    output_clusters += cluster_dict.values()

    output_cluster_id = 0
    output_cluster_ids = [0] * (x_result + 1)
    for cluster in sorted(output_clusters):
        for i in cluster:
            output_cluster_ids[i] = output_cluster_id
        output_cluster_id += 1

    return output_cluster_ids



if __name__ == '__main__':
    path = 'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_Irinotecan_10mg.txt'
    #path = 'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_Oxaliplatin_10mg.txt'
    drug = path.split("th06_")[1].split(".")[0]
    print(f'# drug name: {drug}')
    data_ECv = load_data(path)

    input_data = reshape_inputmatrix(data_ECv)

    result, labels = clustering(input_data)


    savepath = f'resultD_RioEJCA2017/Dendrogram/Dendrogram_RioEJCA2017_SelectedSamples_{drug}.png'
    plot(result, labels, drug, savepath)


    savepath2 = f'resultD_RioEJCA2017/Threshold/Threshold_RioEJCA2017_SelectedSamples_{drug}.png'
    draw_threshold_dependency(result, savepath2)

    # get cluster num & matrix
    clusterIDs = get_cluster_by_number(result, 10)
    print(clusterIDs)


# %%
