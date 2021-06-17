# Author: S.Inoue
# Date: 6/16/2021
# Updated: 6/16/2021
# Project: CTOS folfoliox project
# dataset: RioEJCA2017
# Script: clustering for SurvivalAnalysis

#%%
# import module
from scipy.cluster.hierarchy import linkage, dendrogram
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter


#%%
def load_data(file):
    data = pd.read_table(file, sep='\t', header=0, index_col=0)
    print(f'[LOAD] {os.path.basename(file)}')
    print(f'# input matrix: {data.shape}')
    return data

#%%
def reshape_inputmatrix(data_ClusterIDs_, data_survival):
    data_ClusterIDs_['sampleID'] = [list(data_ClusterIDs_.index)[i].split(':')[0] for i in range(len(data_ClusterIDs_))]
    data_ClusterIDs = data_ClusterIDs_[["clusterID", "sampleID"]]
    reshaped_data = pd.merge(data_ClusterIDs, data_survival,left_on='sampleID',right_index=True).drop("regimen", axis=1)
    return reshaped_data


#%%
def data_visualization(data_, savepath, drug):
    data = data_.sort_values(by="os", ascending=False)
    n_patients = data.shape[0]
    patients = np.arange(n_patients)
    fig, ax = plt.subplots(figsize=(8, 6))
    blue, _, red = sns.color_palette()[:3]
    ax.hlines(patients[data["os censored"].values == 1], 0, data[data["os censored"].values == 1]["os"], color='red', label='Death')
    ax.hlines(patients[data["os censored"].values == 0], 0, data[data["os censored"].values == 0]["os"], color='blue', label='Recovered')
    #ax.scatter(data[data["os censored"].values == 1]['os'], patients[data['os censored'].values == 1], color='k', zorder=10, label='Censored')
    ax.set_title(f"Days since hospitalization [data = RioEJCA2017, drug = {drug}")
    ax.set_xlabel('Days since hospitalization')
    ax.set_ylabel('Case No.')
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1.0))
    plt.savefig(savepath, dpi=100, format='png', bbox_inches="tight")  # save
    print(f'[SAVE]: {savepath}')
    return

#%%
def survival(data, savepath,drug, k):
    df_km = data

    # Create a kmf object
    kmf = KaplanMeierFitter()

    # Fit the data into the model
    #ax = plt.subplot(111)
    fig, ax = plt.subplots(figsize=(8, 6))
    ## cal by each cluster
    for l in range(data["clusterID"].nunique()):
        durations = df_km[df_km["clusterID"] == l]["os"].values
        event_observed = df_km[df_km["clusterID"] == l]["os censored"].values
        kmf.fit(durations, event_observed, label=f'cluster {l+1}')
        kmf.plot(ax=ax)
        #kmf.plot(at_risk_counts=True)
    # label
    ax.set_xlabel("Timeline (days)")
    ax.set_ylabel("Probability event (Recovered or Death) has not occurred")
    ax.set_title(f"KaplanMeier: Survival Analysis [data = RioEJCA2017, drug = {drug}, k = {k}")
    plt.savefig(savepath, dpi=100, format='png', bbox_inches="tight")  # save
    print(f'[SAVE]: {savepath}')
    return


if __name__ == '__main__':
    # load data
    #path = 'data_RioEJCA2017/ClusterIDs_RioEJCA2017_Irinotecan_10mg_k=2.txt'
    #path = 'data_RioEJCA2017/ClusterIDs_RioEJCA2017_Oxaliplatin_10mg_k=2.txt'
    #path = 'data_RioEJCA2017/ClusterIDs_RioEJCA2017_Irinotecan_10mg_k=3.txt'
    path = 'data_RioEJCA2017/ClusterIDs_RioEJCA2017_Oxaliplatin_10mg_k=3.txt'
    drug = path.split("2017_")[1].split("_k=")[0]
    print(f'# drug name: {drug}')
    data_ClusterIDs_ = load_data(path)

    path2 = "data_RioEJCA2017/pfsos_RioEJCA2017.txt"
    data_survival = load_data(path2)

    #
    input_data = reshape_inputmatrix(data_ClusterIDs_, data_survival)

    #
    savepath = f"resultD_RioEJCA2017/SurvivalAnalysis/VisualiveOS_RioEJCA2017_{drug}.png"
    data_visualization(input_data, savepath, drug)

    #
    k = input_data["clusterID"].nunique()
    savepath = f"resultD_RioEJCA2017/SurvivalAnalysis/SurvivalAnalysis_RioEJCA2017_{drug}_k={k}.png"
    survival(input_data, savepath, drug, k)
