# Author: S.Inoue
# Date: 12/21/2020
# Updated: 12/22/2020
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
    # [LOAD] ECv data
    file_1 = 'data/ECv_network_result_0.15.txt'
    data_ecv_ = pd.read_table(file_1, sep='\t', header=0)
    print(f'[LOAD]: {file_1}, input matrix: {data_ecv_.shape}')

    # [LOAD] tumor growth rate data
    file_2 = 'result/txt/IDEA1_1/tumor_growth_rate.pickle'
    with open(file_2, 'rb') as f:
        list_tgr_ = pickle.load(f)
        print(f'[LOAD]: {file_2}, list length: {len(list_tgr_)}')

    # [LOAD] drug name list
    file_3 = 'result/txt/drug_index.csv'
    data_drug_name = pd.read_csv(file_3, sep=',', header=0, index_col=0)
    print(f'[LOAD]: {file_3}, input matrix: {data_drug_name.shape}')

    return data_ecv_, list_tgr_, data_drug_name


# reshape
def reshape(data_ecv_, list_tgr_):
    # reshape data_ecv
    data_ecv = data_ecv_.loc[:, ['Parent',
                                 'Child',
                                 'ECv:C97-float:8',  # CTOS_line = 1
                                 'ECv:C166-float:21',  # CTOS_line = 2
                                 'ECv:C86-float:17',  # CTOS_line = 3
                                 'ECv:C111-foat:18',  # CTOS_line = 4
                                 'ECv:C45-float:5',  # CTOS_line = 5
                                 'ECv:C48-float:6',  # CTOS_line = 6
                                 'ECv:C138-float:20',  # CTOS_line = 7
                                 'ECv:CB3-float:22',  # CTOS_line = 8
                                 'ECv:C75-float:7',  # CTOS_line = 9
                                 'ECv:C132-float:19'  # CTOS_line = 10
                                 ]]  # 必要なedgeを抽出

    # reshape list_tgr
    list_tgr_d28 = []
    for i in range(len(list_tgr_)):
        tgr_d28_ = list_tgr_[i].iloc[:, 8:]
        list_tgr_d28.append(tgr_d28_)

    print('[INFO] reshape the data')
    return data_ecv, list_tgr_d28


# calculate correlation between ECv and TGR, and plot distribution
def corrplot(data_ecv, data_drug_name, list_tgr_d28):
    list_corr_ = []
    for d, i in itertools.product(range(len(data_drug_name['drug_name'])), range(len(data_ecv))):
        # array tumor growth rate
        array_tgr = np.array(list_tgr_d28[d], dtype='float64').reshape(
            10,)  # drug 0~6 → 1
        # array ecv
        array_ecv = np.array(
            data_ecv.iloc[i, 2:], dtype='float64')  # ecv 0~78858 → 1
        # corr
        corr = np.corrcoef(array_tgr, array_ecv)[0, 1]
        # append
        list_corr_.append(corr)
    #reshape
    list_corr = np.array(list_corr_).reshape(len(data_drug_name['drug_name']), -1)

    for l in range(len(list_corr)):
        # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
        colorlist = ["#999999", "#9BC99B", "#3FBF3F",
                    "#9BA7C9", "#3F5FBF", "#D88C8C", "#BF3F3F"]
        # plot
        plt.figure(figsize=(30, 15))
        plt.rcParams["font.family"] = 'sans-serif'
        plt.rcParams['font.size'] = 25
        plt.title(
            f'Distribution of Correlation Coefficients [{data_drug_name["drug_name"][l]}]', fontsize=30)
        plt.xlabel(
            'Correlation coefficients between tumor growth rate [day28] and every ECv', fontsize=25)
        plt.ylabel('Count', fontsize=25)
        sns.distplot(list_corr[l], color=colorlist[l], kde=False)
        # save
        savepath = f'./result/fig/IDEA1_2/1_corrplot/fig_corrplot_{data_drug_name["drug_name"][l]}.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')
    return list_corr


def scatterplot_maxcorrpair(data_drug_name, list_corr, data_ecv, list_tgr_d28):
    for d in range(len(data_drug_name['drug_name'])):
        # information
        ecv_index = list(list_corr[d]).index(max(list_corr[d]))
        drug = data_drug_name["drug_name"][d]
        parent = data_ecv.iloc[ecv_index, :2][0]
        child = data_ecv.iloc[ecv_index, :2][1]
        print(
            f'drug = {drug}, index = {ecv_index}, parent = {parent}, child = {child}')

        # array tumor growth rate & ecv
        array_tgr = np.array(list_tgr_d28[d], dtype='float64').reshape(10,)
        array_ecv = np.array(data_ecv.iloc[ecv_index, 2:], dtype='float64')
        # corr
        corr = np.corrcoef(array_tgr, array_ecv)[0, 1]

        # scatter plot
        plt.figure(figsize=(30, 10))
        plt.rcParams["font.family"] = 'sans-serif'
        plt.rcParams['font.size'] = 25
        plt.title(
            f'scatterplot [{drug}, parent = {parent}, child = {child}]', fontsize=30)
        plt.xlim([0, 16])
        #plt.ylim([0,16])
        plt.xlabel('tumor growth rate [day28]', fontsize=25)
        plt.ylabel('ECv', fontsize=25)
        plt.text(14.2, min(array_ecv),
                 f"corr : {round(max(list_corr[d]), 2)}", size=30)

        data_name = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
        colorlist = ["#999999", "#9BC99B", "#3FBF3F",
                     "#9BA7C9", "#3F5FBF", "#D88C8C", "#BF3F3F"]

        for (x, y, k) in zip(array_tgr, array_ecv, data_name):
            plt.scatter(x, y, color=colorlist[d], s=200)
            plt.annotate(k, xy=(x+0.1, y), size=25)

        # save
        savepath = f'./result/fig/IDEA1_2/2_scatterplot_maxcorrpair/fig_scatterplot_{data_drug_name["drug_name"][d]}.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')
    return

# calculate range between max-min value of ecv
def cal_ecv_maxmin(data_ecv):
    list_ecv_maxmin = []
    for i in range(len(data_ecv)):
        ecv_max = data_ecv.iloc[i, 2:].max()
        ecv_min = data_ecv.iloc[i, 2:].min()
        ecv_maxmin = abs(ecv_max - ecv_min)
        list_ecv_maxmin.append(ecv_maxmin)
    # figure
    plt.figure(figsize=(30, 15))
    plt.rcParams["font.family"] = 'sans-serif'
    plt.rcParams['font.size'] = 25
    plt.title(f'Distribution of range between max-min value of ECv', fontsize=30)
    plt.xlabel('range between max-min value of ECv', fontsize=25)
    plt.ylabel('Count', fontsize=25)
    sns.distplot(list_ecv_maxmin, color="#9BC99B", kde=False)
    # save
    savepath = f'./result/fig/IDEA1_2/3_scatterplot_range_ecv/fig_scatterplot_range_ecv.png'
    plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
    print(f'[SAVE]: {savepath}')
    return list_ecv_maxmin

# make matrix [parent, child, corr, range_ecv]
def corr_matrix(list_corr, list_ecv_maxmin, data_ecv, data_drug_name):
    for i in range(len(list_corr)):
        data = data_ecv.iloc[:, :2]
        data['CorrelationCoefficients'] = list_corr[i]
        data['range_ecv'] = list_ecv_maxmin
        # save
        savepath = f'./result/txt/IDEA1_2/1_CorrelationCoefficients/CorrelationCoefficients_{data_drug_name["drug_name"][i]}.txt'
        data.to_csv(savepath, sep='\t', index=False)
    return



if __name__ == '__main__':
    # load the data
    data_ecv_, list_tgr_, data_drug_name = data_load()
    # reshape the data
    data_ecv, list_tgr_d28 = reshape(data_ecv_, list_tgr_)
    # corr plot
    list_corr = corrplot(data_ecv, data_drug_name, list_tgr_d28)
    # scatter plot
    scatterplot_maxcorrpair(data_drug_name, list_corr, data_ecv, list_tgr_d28)
    # calculate range ecv
    list_ecv_maxmin = cal_ecv_maxmin(data_ecv)
    # corr matrix
    corr_matrix(list_corr, list_ecv_maxmin, data_ecv, data_drug_name)
