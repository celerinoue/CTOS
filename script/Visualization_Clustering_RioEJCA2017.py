# Author: S.Inoue
# Date: 6/14/2021
# Updated: 6/14/2021
# Project: CTOS folfoliox project
# dataset: RioEJCA2017
# Script: clustering for lifelines

#%%
# import module
from scipy.cluster.hierarchy import linkage, dendrogram
import numpy as np
import pandas as pd
import os
import seaborn as sns
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
def clustering(data, data_pfsos):
    labels = list(data.index)
    result = linkage(data.iloc[:, :],
                    #metric = 'braycurtis',
                    #metric = 'canberra',
                    #metric = 'chebyshev',
                    #metric = 'cityblock',
                    #metric='correlation',
                    #metric = 'cosine',
                    metric = 'euclidean',
                    #metric = 'hamming',
                    #metric = 'jaccard',
                    #method= 'single')
                    #method='average')
                    method = 'ward')
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


def get_cluster_by_number(result, number, data, savepath):
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

    # make matrix
    data["clusterID"] = output_cluster_ids
    # savelist
    data.to_csv(savepath, sep='\t', index=True)
    print(f'[SAVE]: {savepath}')

    return output_cluster_ids


def heatmap(data, data_pfsos, labels, savepath, drug):
    # set res
    list_res_ = [labels[i].split(":")[0] for i in range(len(labels))]
    list_res = []
    for k in range(len(list_res_)):
        res = data_pfsos[data_pfsos.index ==list_res_[k]]["response status"].to_list()
        list_res.append(res)
    list_res = list(itertools.chain.from_iterable(list_res))
    # set color
    list_colors = []
    for i in list_res:
        l = i.replace("NR", "blue").replace("R", "red")
        list_colors.append(l)

    # heatmap
    sns.clustermap(data, method='ward', metric='euclidean', row_cluster=True, col_cluster=True, row_colors=list_colors)
    plt.title(f"heatmap [drug = {drug}]")
    plt.savefig(savepath, dpi=100, format='png', bbox_inches="tight")  # save
    print(f'[SAVE]: {savepath}')
    return


if __name__ == '__main__':
    # pfsos
    path2 = "data_RioEJCA2017/pfsos_RioEJCA2017.txt"
    data_pfsos = load_data(path2)
    # ECv
    #path = 'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_Irinotecan_10mg_v2.txt'
    #path = 'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_Oxaliplatin_10mg_v2.txt'
    path = 'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_Irinotecan_10mg_v3.txt'
    #path = 'data_2ndFeatureExtractedDataset/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06/2ndFeatureExtractedDataset_RioEJCA2017_ECv_th06_Oxaliplatin_10mg_v3.txt'

    drug = path.split("th06_")[1].split("_v")[0]
    print(f'# drug name: {drug}')
    data_ECv = load_data(path)

    input_data = reshape_inputmatrix(data_ECv)

    result, labels = clustering(input_data, data_pfsos)



    savepath = f'resultD_RioEJCA2017/Dendrogram/Dendrogram_RioEJCA2017_SelectedSamples_{drug}_v3.png'
    plot(result, labels, drug, savepath)


    savepath2 = f'resultD_RioEJCA2017/Threshold/Threshold_RioEJCA2017_SelectedSamples_{drug}_v3.png'
    draw_threshold_dependency(result, savepath2)

    # get cluster num & matrix
    k = 3
    savepath3 = f'data_RioEJCA2017/ClusterIDs_RioEJCA2017_{drug}_k={k}_v3.txt'
    clusterIDs = get_cluster_by_number(result, k, input_data, savepath3)
    print(clusterIDs)

    # heatmap
    savepath4 = f'resultD_RioEJCA2017/Heatmap/Heatmap_RioEJCA2017_SelectedSamples_{drug}_v3.png'
    heatmap(input_data, data_pfsos, labels, savepath4, drug)
