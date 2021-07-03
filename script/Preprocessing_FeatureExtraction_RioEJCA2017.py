# Author: S. Inoue
# Date: 5/28/2021
# Updated: 05/28/2021
# Project: CTOS
# Scropt: To generate dataset from raw microarray&clinical data for BN input
# Array dataset: RioEJCA2017 (GSE147571) for Survival Analysis

#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# %%

#%%
# data load
def dataload(file_1, file_2):
    # gene name list
    with open(file_1) as f:
        # gene list
        lines = f.read().splitlines()[270:54946]
        gene_list_ = pd.DataFrame([s.split("\t") for s in lines])
        gene_list_.columns = gene_list_.iloc[0]
        gene_list_ = gene_list_.drop(0)
        gene_list = gene_list_[["ID", "Representative Public ID", "Gene Symbol"]]

    # gene expression data
    with open(file_2) as f:
        # data INFO
        lines = f.read().splitlines()[0:28]
        info = pd.DataFrame([s.split("\t") for s in lines])
        info[0] = info[0].str.lstrip('!')
        info[1] = info[1].str.lstrip('\"').str.rstrip('\"')
    with open(file_2) as f:
        # series_matrix_info
        lines_1 = f.read().splitlines()[29:74]
        data_1 = pd.DataFrame([s.split("\t") for s in lines_1])
        data_1[0] = data_1[0].str.lstrip('!')
        for i in range(len(data_1.columns)):
            data_1[i] = data_1[i].str.lstrip('\"').str.rstrip('\"')
        data_1.columns = data_1.iloc[1]
        #data_1 = data_1.drop()
    with open(file_2) as f:
        # series_matrix_table
        lines_2 = f.read().splitlines()[75:54676]
        data_2 = pd.DataFrame([s.split("\t") for s in lines_2])
        data_2[0] = data_2[0].str.lstrip('!')
        for i in range(len(data_2.columns)):
            data_2[i] = data_2[i].str.lstrip('\"').str.rstrip('\"')
        data_2.columns = data_2.iloc[0]
        data_2 = data_2.drop(0)

    # reshape =============
    # data_2
    # gene_expression matrix (index = GeneName , allow Duplicate)
    gene_exp__ = pd.merge(data_2, gene_list, left_on='ID_REF', right_on='ID').drop(columns=[
        "ID_REF", "ID", "Representative Public ID"]).rename(columns={'Gene Symbol': 'GeneName'})
    gene_exp_ = gene_exp__.drop("GeneName", axis=1).astype(float)  # 数値変換
    gene_exp_["GeneName"] = gene_exp__["GeneName"]
    gene_exp = gene_exp_.groupby(by="GeneName").mean().drop(index=['']).dropna(how="any")  # 重複削除
    # log2
    #gene_exp = np.log2(gene_exp + 1)[1:]


    # data_1
    char_list = ["center", "sex", "age", "who performance status", "tumor location", "synchronous metastase",
                 "pn", "pt", "regimen", "response category", "response status", "pfs censored", "pfs", "os censored", "os"]
    data_1["Sample_geo_accession"][9:24] = char_list
    data_1 = data_1.set_index("Sample_geo_accession").T
    pfs_os = data_1[["regimen", "os censored", "os", "response status"]]
    pfs_os["regimen"] = [pfs_os["regimen"][i].split(":")[1].replace(' ', '') for i in range(len(pfs_os))]
    pfs_os["os censored"] = [pfs_os["os censored"][i].split(":")[1].replace(' ', '') for i in range(len(pfs_os))]
    #pfs_os["pfs"] = pfs_os["pfs"].astype(float)
    pfs_os["os"] = [pfs_os["os"][i].split(":")[1] for i in range(len(pfs_os))]
    pfs_os["os"] = pfs_os["os"].astype(float)
    # "response status"
    pfs_os["response status"] = [pfs_os["response status"][i].split(":")[1].replace(' ', '') for i in range(len(pfs_os))]
    pfs_os = pfs_os.reset_index().rename(columns={1: 'Sample_name'})

    return gene_exp, pfs_os




def histogram(matrix, filename):
    ## plot histogram per genes ##
    plt.hist(matrix.mean(axis='columns'), alpha=0.4, bins=30, color='lime')
    plt.xlabel('expression')
    plt.ylabel('frequency')
    plt.savefig(filename, dpi=300, format='png')
    plt.clf()
    print(f'[SAVE]: {filename}')




def violinplot(matrix, filename):
    ## violinplot per sample to see data distribution ##
    plt.figure(figsize=(15, 8))
    sns.violinplot(data=matrix, inner='box')
    plt.ylabel('expression')
    plt.xticks(fontsize=9, weight='bold', rotation=70)
    plt.savefig(filename, dpi=300, format='png')
    plt.clf()
    print(f'[SAVE]: {filename}')


    fig = plt.figure(figsize=(8, 6), tight_layout=True)
    ax=fig.add_subplot(1, 1, 1)
    ax.spines['top'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_linewidth(1)




if __name__ == '__main__':
    # load data
    file_1 = 'data/RioEJCA2017/GSE72970_family.soft'
    file_2 = 'data/RioEJCA2017/GSE72970_series_matrix.txt'

    gene_exp, pfs_os = dataload(file_1, file_2)

    print('[INFO] feature extraction completed')

    # save file
    savepath1 = 'data_RioEJCA2017/FeatureExtractedMatrix_RioEJCA2017_v2.txt'
    gene_exp.to_csv(savepath1, mode='w', sep="\t")
    print(f"[SAVE] {savepath1}")

    savepath2 = 'data_RioEJCA2017/pfsos_RioEJCA2017.txt'
    pfs_os.to_csv(savepath2, mode='w', sep="\t", index=False)
    print(f"[SAVE] {savepath2}")


# %%

# %%
