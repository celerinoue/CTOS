# Author: yoshi
# Date: 9/18/2020
# Updated: 11/12/2020
# Project: CTOS folfoliox project
# Script: To perform basic analyses

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from pandas import plotting

def dataload(filename):
    ## load data ##
    data = pd.read_table(filename, sep='\t', header=0, index_col=0)
    print(f'[LOAD]: {filename}')
    print(f'input matrix: {data.shape}')
    return data


def scatterplot(matrix, filename):
    ## Depict scatter plot across all pairwise samples ##
    plt.figure(figsize=(25, 25))
    plotting.scatter_matrix(matrix, diagonal='hist', alpha=0.2, s=1, c='lime')
    plt.savefig(filename, dpi=300, format='png')
    plt.clf()
    print(f'[SAVE]: {filename}')


def histogram(matrix, filename):
    ## plot histogram per genes ##
    plt.hist(matrix.mean(axis='columns'), alpha=0.4, bins=30, color='lime')
    plt.xlabel('expression')
    plt.ylabel('frequency')
    plt.savefig(filename, dpi=300, format='png')
    plt.clf()
    print(f'[SAVE]: {filename}')


def pca(matrix, filename):

    samplenames = list(matrix.columns)

    ## Standard scaling for PCA input (mean:0, variance:1)
    scaler = StandardScaler()
    scaler.fit(matrix)
    data_std = scaler.transform(matrix) #ndarray
    #data_std_df=pd.DataFrame(data_std, columns=sample_names) #dataframe

    ## PCA
    n_components = 2
    pca = PCA(n_components)
    pca.fit(data_std.T)
    pca_output = pca.transform(data_std.T)
    embed = pd.DataFrame(pca_output, columns=['PC1', 'PC2'], index=samplenames)

    plt.scatter(embed['PC1'], embed['PC2'], alpha=0.6, c='mediumblue', s=23)
    for x, y, name in zip(embed['PC1'], embed['PC2'], samplenames):
        plt.text(x, y, name, alpha=0.6, fontsize=7, weight='bold')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig(filename, dpi=300, format='png')
    plt.clf()
    print(f'[SAVE]: {filename}')


def violinplot(matrix, filename):
    ## violinplot per sample to see data distribution ##
    plt.figure(figsize=(15, 8))
    sns.violinplot(data=matrix, inner='box')
    plt.ylabel('expression')
    plt.xticks(fontsize=9, weight='bold',rotation=70)
    plt.savefig(filename, dpi=300, format='png')
    plt.clf()
    print(f'[SAVE]: {filename}')


if __name__ == '__main__':

    SCNECdata = dataload('../data/SCNEC_dataset.txt') #theme color: deeppink
    CRCdata = dataload('../data/CRC_dataset.txt') #theme color: lime
    MIXdata = pd.concat([SCNECdata, CRCdata], join='inner', axis=1)
    print(f'MIXdata shape: {MIXdata.shape}')

    ## CRC: select float condition (theme color: mediumblue)
    CRCdata_float = CRCdata.loc[:,['C45-float',
                                   'C48-float',
                                   'C75-float',
                                   'C97-float',
                                   'C86-float',
                                   'C111-foat',
                                   'C132-float',
                                   'C138-float',
                                   'C166-float',
                                   'CB3-float',
                                   'HT29',
                                   'SW620']]
    print(f'CRC-float shape: {CRCdata_float.shape}')

    ## CRC: select non condition
    CRCdata_non = CRCdata.loc[:,['C45',
                                 'C48',
                                 'C75',
                                 'C97',
                                 'C86',
                                 'C111',
                                 'C132',
                                 'C138',
                                 'C166',
                                 'CB3',
                                 'HT29',
                                 'SW620']]
    print(f'CRC-non shape: {CRCdata_non.shape}')

    scatterplot(SCNECdata, '../data/SCNEC_scatterplot.png')
    scatterplot(CRCdata, '../data/CRC_scatterplot.png')

    histogram(SCNECdata, '../data/SCNEC_histogram.png')
    histogram(CRCdata, '../data/CRC_histogram.png')

    pca(SCNECdata, '../data/SCNEC_pca.png')
    pca(CRCdata, '../data/CRC_pca.png')
    pca(MIXdata, '../data/MIXdata_pca.png') 
    pca(CRCdata_float, '../data/CRC_float_pca.png')
    pca(CRCdata_non, '../data/CRC_non_pca.png')

    violinplot(SCNECdata, '../data/SCNEC_violinplot.png')
    violinplot(CRCdata, '../data/CRC_violinplot.png')
