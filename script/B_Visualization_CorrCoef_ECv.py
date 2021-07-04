# Author: S.Inoue
# Date: 12/21/2020
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
    # [LOAD] ECv data
    file_1 = 'BayesianNetworkEstimation/CRC/ECv_network_result_0.15.txt'
    data_ecv_ = pd.read_table(file_1, sep='\t', header=0)
    print(f'[LOAD]: {file_1}, input matrix: {data_ecv_.shape}')

    # [LOAD] tumor growth rate data (TGR data, drug1~7)
    file_2 = 'data/data_TumorGrowthRate/data_TumorGrowthRate.pickle'
    with open(file_2, 'rb') as f:
        list_tgr_ = pickle.load(f)
        print(f'[LOAD]: {file_2}, list length: {len(list_tgr_)}')

    # [LOAD] drug name list
    file_3 = 'data/data_TumorGrowthRate/drug_index.csv'
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


# ============================== calculate ==============================

# calculate correlation between [ECv & TGR]
def calculate_corr_ecv_tgr(data_ecv, data_drug_name, list_tgr_d28):
    list_corr_ = []
    for d, i in itertools.product(range(len(data_drug_name['drug_name'])), range(len(data_ecv))):
        # array TGR
        array_tgr = np.array(list_tgr_d28[d], dtype='float64').reshape(10,)  # drug 0~6 → 1
        # array gene
        array_ecv = np.array(data_ecv.iloc[i, 2:], dtype='float64')  # ecv 0~78858 → 1
        # corr
        corr = np.corrcoef(array_tgr, array_ecv)[0, 1]
        # append
        list_corr_.append(corr)
    #reshape
    list_corr = np.array(list_corr_).reshape(len(data_drug_name['drug_name']), -1)

    return list_corr

# calculate range of ECv max-min value and save matrix
def cal_ecv_range(data_ecv, data_drug_name, list_corr):
    # calculate range of ECv max-min value
    list_ecv_range = []
    for i in range(len(data_ecv)):
        ecv_max = data_ecv.iloc[i, 2:].max()  # LINE1~10の中でのECvの最大値
        ecv_min = data_ecv.iloc[i, 2:].min()  # LINE1~10の中でのECvの最小値
        ecv_range = abs(ecv_max - ecv_min)  # 絶対値(最大値 - 最小値)
        list_ecv_range.append(ecv_range)
    print(f'[INFO] make list of ECv max-min range [shape = {len(list_ecv_range)}]')

    # make matrix [1.parent, 2.child, 3.corr, 4.range_ecv]
    for d in range(len(data_drug_name)):
        matrix = data_ecv.iloc[:, :2]  # 1,2列目
        matrix['CorrelationCoefficients'] = list_corr[d]
        matrix['range_ecv'] = list_ecv_range  # 4列目
        # save
        savepath = f'data/data_CorrCoef/ECv/CorrelationCoefficients_{data_drug_name["drug_name"][d]}.txt'
        matrix.to_csv(savepath, sep='\t', index=False)

    return list_ecv_range


# ============================== figure ==============================

# plot distribution of correlation values[ECv & TGR]
def plot_distribution_corr_ecv_tgr(list_corr, data_drug_name):
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
        plt.title(f'Distribution of Correlation Coefficients between ECv and Tumor Growth Rate(day28) [drug = {data_drug_name["drug_name"][l]}]', fontsize=30)
        plt.xlabel('Correlation coefficients between ECv and TGR(day28)', fontsize=25)
        plt.ylabel('Frequency', fontsize=25)
        sns.distplot(list_corr[l], color=colorlist[l], kde=False)
        # save
        savepath = f'resultB_CorrCoef/range_ECv/dist_CorrCoef_ECv/fig_corrplot_{data_drug_name["drug_name"][l]}.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')

    return


# Plot a scatter plot of pairs with the largest correlation coefficient
def plot_scatter_maxcorr(data_ecv, data_drug_name, list_corr):
    # 各薬剤ごとに
    for d in range(len(data_drug_name)):
        # 相関係数最大のecvのindexと値, parent, childの名前を取得=========
        maxcorr_index = list(list_corr[d]).index(max(list_corr[d]))  # indexを取得
        drug = data_drug_name["drug_name"][d]
        parent = data_ecv.iloc[maxcorr_index, :2][0]  # parentの名前を取得
        child = data_ecv.iloc[maxcorr_index, :2][1]  # childの名前を取得
        tgr_value = np.array(list_tgr_d28[d], dtype='float64').reshape(
            10,)  # CTOS1~10のtgrの値
        ecv_value = np.array(
            data_ecv.iloc[maxcorr_index, 2:], dtype='float64')  # CTOS1~10のecvの値
        print(
            f'drug = {drug}, index = {maxcorr_index}, parent = {parent}, child = {child}')

        # scatter plot==================
        plt.figure(figsize=(30, 10))
        plt.rcParams["font.family"] = 'sans-serif'
        plt.rcParams['font.size'] = 25
        plt.title(
            f'scatterplot [{drug}, parent = {parent}, child = {child}]', fontsize=30)
        plt.xlim([0, 16])
        #plt.ylim([0,16])
        plt.xlabel('tumor growth rate [day28]', fontsize=25)
        plt.ylabel('ECv', fontsize=25)
        plt.text(14.2, min(ecv_value),
                 f"corr : {round(max(list_corr[d]), 2)}", size=30)

        CTOS_LINE = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
        colorlist = ["#999999",  # gray
                     "#9BC99B", "#3FBF3F",  # green
                     "#9BA7C9", "#3F5FBF",  # blue
                     "#D88C8C", "#BF3F3F"  # red
                     ]

        for (x, y, k) in zip(tgr_value, ecv_value, CTOS_LINE):
            plt.scatter(x, y, color=colorlist[d], s=200)
            plt.annotate(k, xy=(x+0.1, y), size=25)

        # save figure===============
        savepath = f'resultB_CorrCoef/range_ECv/scat_maxCorrCoef_ECv/fig_scatterplot_{data_drug_name["drug_name"][d]}.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')

    return


# plot distribution of range [ECv & TGR]
def plot_distribution_ecv_range(list_ecv_range):
    # figure
    plt.figure(figsize=(30, 15))
    plt.rcParams["font.family"] = 'sans-serif'
    plt.rcParams['font.size'] = 25
    plt.title(f'Distribution of range between max-min value of ECv', fontsize=30)
    plt.xlabel('range between max-min value of ECv', fontsize=25)
    plt.ylabel('Frequency', fontsize=25)
    sns.distplot(list_ecv_range, color="#9BC99B", kde=False)
    # save
    savepath = f'resultB_CorrCoef/range_ECv/dist_CorrCoef_ECv/fig_distribution_range_ecv.png'
    plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
    print(f'[SAVE]: {savepath}')
    return


if __name__ == '__main__':
    # load data
    data_ecv_, list_tgr_, data_drug_name = data_load()
    # reshape the data
    data_ecv, list_tgr_d28 = reshape(data_ecv_, list_tgr_)

    # calculate correlation between [ECv & TGR]
    list_corr = calculate_corr_ecv_tgr(data_ecv, data_drug_name, list_tgr_d28)
    # calculate range of ECv max-min value and save matrix
    list_ecv_range = cal_ecv_range(data_ecv, data_drug_name, list_corr)

    # plot distribution of correlation values[ECv & TGR]
    plot_distribution_corr_ecv_tgr(list_corr, data_drug_name)
    # plot a scatter plot of pairs with the largest correlation coefficient
    plot_scatter_maxcorr(data_ecv, data_drug_name, list_corr)
    # plot distribution of range [ECv & TGR]
    plot_distribution_ecv_range(list_ecv_range)
