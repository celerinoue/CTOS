# Author: S.Inoue
# Date: 01/21/2021
# Updated: 03/18/2021
# Project: CTOS folfoliox project
# Script: To perform basic analyses

# import module
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

# load the data
def data_load():
    # [LOAD] gene exprssion data
    file_1 = 'BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt'
    data_gene_ = pd.read_table(file_1, sep='\t', header=0)
    print(f'[LOAD]: {file_1}, input matrix: {data_gene_.shape}')

    # [LOAD] tumor growth rate data (TGR data, drug1~7)
    file_2 = 'data_TGR/TGR.pickle'
    with open(file_2, 'rb') as f:
        list_tgr_ = pickle.load(f)
        print(f'[LOAD]: {file_2}, list length: {len(list_tgr_)}')

    # [LOAD] drug name list
    file_3 = 'data_TGR/drug_index.csv'
    data_drug_name = pd.read_csv(file_3, sep=',', header=0, index_col=0)
    print(f'[LOAD]: {file_3}, input matrix: {data_drug_name.shape}')

    return data_gene_, list_tgr_, data_drug_name


# reshape the data
def reshape(data_gene_, list_tgr_):
    # reshape data_gene 必要なcolumnのみを抽出
    data_gene = data_gene_.loc[:, ['GeneName',
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

    # reshape list_tgr 28日目のTGRを取得
    list_tgr_d28 = []
    for i in range(len(list_tgr_)):
        tgr_d28_ = list_tgr_[i].iloc[:, 8:]
        list_tgr_d28.append(tgr_d28_)

    print('[INFO] reshape the data')
    return data_gene, list_tgr_d28


# calculate correlation between [gene & TGR]
def calculate_corr_gene_tgr(data_gene, data_drug_name, list_tgr_d28):
    list_corr_ = []
    for d, i in itertools.product(range(len(data_drug_name['drug_name'])), range(len(data_gene))):
        # array TGR
        array_tgr = np.array(list_tgr_d28[d], dtype='float64').reshape(10,)  # drug 0~6 → 1
        # array gene
        array_gene = np.array(data_gene.iloc[i, 1:], dtype='float64')  # gene 0~78858 → 1
        # corr
        corr = np.corrcoef(array_tgr, array_gene)[0, 1]
        # append
        list_corr_.append(corr)
    #reshape
    list_corr = np.array(list_corr_).reshape(len(data_drug_name['drug_name']), -1)

    return list_corr


# calculate range of gene expression max-min value and save matrix
def cal_gene_range(data_gene, data_drug_name, list_corr):
    # calculate range of gene expression max-min value
    list_gene_range = []
    for i in range(len(data_gene)):
        gene_max = data_gene.iloc[i, 1:].max()  # LINE1~10の中でのgeneの最大値
        gene_min = data_gene.iloc[i, 1:].min()  # LINE1~10の中でのgeneの最小値
        gene_range = abs(gene_max - gene_min)  # 絶対値(最大値 - 最小値)
        list_gene_range.append(gene_range)
    print(
        f'[INFO] make list of gene expression max-min range [shape = {len(list_gene_range)}]')

    # make matrix [1.gene, 3.corr, 4.range_gene]
    for d in range(len(data_drug_name)):
        matrix = pd.DataFrame(data_gene.iloc[:, 0])  # 1列目
        matrix['CorrelationCoefficients'] = list_corr.copy()[d]
        matrix['range_gene'] = list_gene_range # 4列目
        # save
        savepath = f'resultA_CorrCoef/range_GeneExp/CorrelationCoefficients_{data_drug_name["drug_name"][d]}.txt'
        matrix.to_csv(savepath, sep='\t', index=False)

    return list_gene_range


# ====================== plot =============================


# plot distribution of correlation values[gene & TGR]
def plot_distribution_corr_gene_tgr(list_corr, data_drug_name):
    for l in range(len(list_corr)):
        # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
        colorlist = ["#999999",  # gray
                     "#9BC99B", "#3FBF3F",  # green
                     "#9BA7C9", "#3F5FBF",  # blue
                     "#D88C8C", "#BF3F3F"  # red
                     ]
        # plot
        plt.figure(figsize=(30, 15))
        plt.rcParams["font.family"] = 'sans-serif'
        plt.rcParams['font.size'] = 25
        plt.title(f'Distribution of Correlation Coefficients between Gene Expression and Tumor Growth Rate(day28) [drug = {data_drug_name["drug_name"][l]}]', fontsize=30)
        plt.xlabel('Correlation coefficients between Gene Expression and TGR(day28)', fontsize=25)
        plt.ylabel('Frequency', fontsize=25)
        sns.distplot(list_corr[l], color=colorlist[l], kde=False)
        # save
        savepath = f'resultA_CorrCoef/range_GeneExp/dist_CorrCoef_GeneExp/fig_corrplot_{data_drug_name["drug_name"][l]}.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')

    return


# plot distribution of range [gene & TGR]
def plot_distribution_gene_range(list_gene_range):
    # figure
    plt.figure(figsize=(30, 15))
    plt.rcParams["font.family"] = 'sans-serif'
    plt.rcParams['font.size'] = 25
    plt.title(f'Distribution of Difference between max-min gene expression values', fontsize=30)
    plt.xlabel('Difference between the maximum and minimum gene expression values', fontsize=25)
    plt.ylabel('Frequency', fontsize=25)
    sns.distplot(list_gene_range, color="#9BC99B", kde=False)
    # save
    savepath = f'resultA_CorrCoef/range_GeneExp/dist_CorrCoef_GeneExp/fig_distribution_range_GeneExp.png'
    plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
    print(f'[SAVE]: {savepath}')

    return


if __name__ == '__main__':
    # load the data
    data_gene_, list_tgr_, data_drug_name = data_load()
    # reshape the data
    data_gene, list_tgr_d28 = reshape(data_gene_, list_tgr_)

    # calculate correlation between [gene & TGR]
    list_corr = calculate_corr_gene_tgr(data_gene, data_drug_name, list_tgr_d28)
    # calculate range of gene expression max-min value and save matrix
    list_gene_range = cal_gene_range(data_gene, data_drug_name, list_corr)

    # plot distribution of correlation values[gene & TGR]
    plot_distribution_corr_gene_tgr(list_corr, data_drug_name)
    # plot distribution of range [gene & TGR]
    plot_distribution_gene_range(list_gene_range)
