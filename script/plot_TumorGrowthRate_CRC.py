# Author: S.Inoue and yoshi
# Date: 11/5/2020
# Updated: 11/16/2020
# Project: CTOS folfoliox
# Script: To plot and analyze in vivo eperiment result

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import os
import sys
import itertools
from scipy import stats


def data_load(filename):
    data = pd.read_csv(filename)
    print(f'[LOAD]: {filename}')
    print(f'input matrix: {data.shape}')
    return data

# Data Transforming
def data_transform(data):
    ## patient_id 4
    # 7.Irinotecan_10mg → Irinotecan_3mg, 8.Irinotecan_30mg → Irinotecan_10mg
    data[data["patient_id"] == 4] = data[data["patient_id"] == 4].replace({7:6, 8:7, 'Irinotecan_10mg':'Irinotecan_3mg', 'Irinotecan_30mg':'Irinotecan_10mg'})

    ## Convert all value to multiplier
    day_list = ["04", "07", "11", "14", "18", "21", "25", "28"] # day list without 00
    for l in day_list:
        data[l] = data[l] / data["00"] # multiplier
    data["00"] = 1  # day 00

    ## Transform the data format
    pat_list = data["patient_id"].drop_duplicates().tolist()  # patient_idの範囲
    drug_list = data["drug_id"].drop_duplicates().tolist()  # drug_idの範囲
    sub_list = data["subject_id"].drop_duplicates().tolist()  # subject_idの範囲
    l_ = []
    for pat, drg, sub in itertools.product(pat_list, drug_list, sub_list):
        data_row = data[(data['patient_id'] == pat) & (data['drug_id'] == drg) & (
        data['subject_id'] == sub)]  # Specify a column from data
        # define (day, value) from data
        value = data_row.iloc[:, 4:].T.reset_index()
        value.columns = ['day', 'value']  # column names of (day, value)
        # get index (patient_idm drug_id, drug_name, subject_id) from data
        index = data_row.iloc[:, :4].reset_index(drop=True)
        matrix = pd.concat([index, value], axis=1).fillna(method='ffill')
        l_.append(matrix)
    data_val = pd.concat(l_, axis=0).reset_index(drop=True)
    data_val = data_val.astype({'patient_id': int, 'drug_id': int, 'subject_id': int, 'day': int}) #intに変換

    print('[INFO] TRANSFORMED THE DATA!!')
    print(f'data shape: {data_val.shape}')
    return data_val, pat_list, drug_list


# write factorplot [mean, median]
def fig_effects_per_patient(data, pat_num, estimator, ci):
    # option
    '''
    Args:
        data :
        pat_num : patient number [1~10]
        estimator : (optional) Statistical function to estimate [mean or median]
        ci : (optional) Size of confidence intervals to draw around estimated values [float or “sd”]
    '''
    #Args [estimator]
    if estimator == 'median':
        e = np.median
    else:
        estimator = 'mean'
        e = np.mean

    #Args [ci]
    if type(ci) is int:
        ci = ci
    else:
        ci = 'sd'

    #Args [data, pat_num]
    if type(pat_num) is int:
        data = data[data["patient_id"] == pat_num]
        pat_num = pat_num
        savepath = f'./result/fig/{estimator}_plot/fig_{estimator}_ci_{ci}_patient_{pat_num}.png'
    else:
        data = data
        pat_num = 'All [1~10]'
        savepath = f'./result/fig/{estimator}_plot/fig_{estimator}_ci_{ci}_patient_all.png'

    # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
    colorlist = ["#999999", "#9BC99B", "#3FBF3F",
                 "#9BA7C9", "#3F5FBF", "#D88C8C", "#BF3F3F"]
    sns.set_palette(colorlist, 7)
    #  Defaults are size=5, aspect=1
    sns.catplot(x='day', y='value', data=data, kind = 'point', hue='drug_name',
                   estimator=e, ci=ci, height=8, aspect=2)
    # supplementary wire [pink]
    plt.hlines(1, -0.3, 8.3, '#D88C8C', linestyles='dashed')
    sns.set_style("darkgrid", {'grid.linestyle': '--'})
    plt.title(f'Patient {pat_num}, {estimator}, CI = {ci}', fontsize=20)

    plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
    print(f'[SAVE]: {savepath}')


def cal_pvalue(data, path):
    # make list
    pat_list = data["patient_id"].drop_duplicates().tolist()  # patient_idの範囲
    drg_list = data["drug_id"].drop_duplicates().tolist()  # drug_idの範囲
    drg_name_list = data["drug_name"].drop_duplicates().tolist()
    day_list = data["day"].drop_duplicates().tolist()
    day_list.remove(0)  # drop '0' day

    # prepare dataset to calculate p-value
    list_t = []
    for pat, day, drg in itertools.product(pat_list, day_list, drg_list):
        v = data[(data['patient_id'] == pat) & (data['day'] == day)
                 & (data['drug_id'] == drg)]['value']
        list_t.append(v)
        a = np.array(list_t)
    data_reshape = a.reshape([80, 7, 6])

    # calculate p-value
    list_p = []
    for m in range(0, 80):
        for l in range(1, 7):  # t-testの回数
            p = stats.ttest_rel(data_reshape[m][0], data_reshape[m][l])[1]
            list_p.append(p)
    data_pval_p = pd.DataFrame(np.reshape(
        list_p, (80, 6)), columns=drg_name_list[1:7])

    # make dataframe drug_name
    data_pval_f = data[["patient_id", "day"]].drop_duplicates()
    data_pval_f = data_pval_f[data_pval_f['day'] != 0].reset_index(drop=True)

    data_pval = pd.concat([data_pval_f, data_pval_p], axis=1)
    data_pval.to_csv(path)
    print(f'[SAVE]: {path}')
    return


def fig_effects_per_drug(data, drg_num, estimator, ci):
    #Args [estimator]
    if estimator == 'median':
        e = np.median
    else:
        estimator = 'mean'
        e = np.mean

    #Args [ci]
    if type(ci) is int:
        ci = ci
    else:
        ci = 'sd'

    #Args [data, drg_num]
    if type(drg_num) is int:
        data = data[data["drug_id"] == drg_num]
        drg_name = data[data["drug_id"] == drg_num]["drug_name"].iloc[1]
        savepath = f'./result/fig/{estimator}_plot/fig_{estimator}_ci_{ci}_drug_{drg_name}.png'
    else:
        data = data
        drg_name = 'All [1~10]'
        savepath = f'./result/fig/{estimator}_plot/fig_{estimator}_ci_{ci}_drug_all.png'

    sns.set_palette('Set3')
    sns.catplot(x='day', y='tumor growth rate', data=data, kind='point',
                hue='patient_id', estimator=e, ci=ci, height=8, aspect=2)
    #  Defaults are size=5, aspect=1
    plt.hlines(1, -0.3, 8.3, '#D88C8C', linestyles='dashed')  # 補助線
    sns.set_style("darkgrid", {'grid.linestyle': '--'})
    plt.title(f'Drug: {drg_name}, {estimator}, CI = {ci}', fontsize=20)

    plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
    print(f'[SAVE]: {savepath}')



if __name__ == '__main__':
    CTOSdata_ = data_load('./data/CTOS_result.csv')  # load the data
    CTOSdata, pat_list, drug_list = data_transform(CTOSdata_)

    # Figure write [per patient]
    for l in pat_list:
        fig_effects_per_patient(CTOSdata, l, 'mean', 'sd')  # mean plot
        fig_effects_per_patient(CTOSdata, l, 'median', 'sd')  # median plot

    for d in drug_list:
        fig_effects_per_drug(CTOSdata, d, 'mean', 'sd')
        fig_effects_per_drug(CTOSdata, d, 'median', 'sd')

    # Figure write [all patient]
    fig_effects_per_patient(CTOSdata, all, 'mean', 'sd')  # mean plot
    fig_effects_per_patient(CTOSdata, all, 'median', 'sd')  # median plot


    # calculate p-value
    cal_pvalue(CTOSdata, "./result/txt/p_value.csv")

