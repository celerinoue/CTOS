# Author: S.Inoue and yoshi
# Date: 11/5/2020
# Updated: 03/18/2021
# Project: CTOS folfoliox
# Script: Labeling Drug Responsiveness of 3 drugs
# Process: Idea1, Step1

#%%
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import itertools
import pickle
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette

#%%
def dataload(filename):
    data = pd.read_csv(filename)
    print(f'[LOAD]: {filename}')
    print(f'input matrix: {data.shape}')
    return data

## Data Transforming ##
def data_transform(data, savepath):
    ## Modify the original data ##
    # patient_id : 7.Irinotecan_10mg → Irinotecan_3mg, 8.Irinotecan_30mg → Irinotecan_10mg
    data[data["patient_id"] == 4] = data[data["patient_id"] == 4].replace(
        {7: 6, 8: 7, 'Irinotecan_10mg': 'Irinotecan_3mg', 'Irinotecan_30mg': 'Irinotecan_10mg'})

    ## Convert all data to ratio ##
    data_ratio = data  # copy
    day_list = ["00", "04", "07", "11", "14", "18", "21", "25", "28"]
    for l in day_list[1:]:  # day list without 00
        data_ratio[l] = data_ratio[l] / data_ratio["00"]  # cal ratio
    data_ratio["00"] = 1  # add day 00 = ratio 1

    ## Transform the data format ##
    patient_num = data_ratio["patient_id"].drop_duplicates().tolist()  # range patient_id
    drug_num = data_ratio["drug_id"].drop_duplicates().tolist()  # range drug_id
    sub_num = data_ratio["subject_id"].drop_duplicates().tolist()  # range subject_id
    l_ = []
    for pat, drg, sub in itertools.product(patient_num, drug_num, sub_num):
        ## define columns from data ##
        data_row = data_ratio[(data_ratio['patient_id'] == pat) & (data_ratio['drug_id'] == drg) & (
            data_ratio['subject_id'] == sub)]
        ## define columns(day, value) from data ##
        value = data_row.iloc[:, 4:].T.reset_index()
        value.columns = ['day', 'value']
        ## get index (patient_idm drug_id, drug_name, subject_id) from data ##
        index = data_row.iloc[:, :4].reset_index(drop=True)
        matrix = pd.concat([index, value], axis=1).fillna(method='ffill')
        l_.append(matrix)
    data_val = pd.concat(l_, axis=0).reset_index(drop=True)
    data_val = data_val.astype({'patient_id': int, 'drug_id': int, 'subject_id': int, 'day': int})  # convert to int

    print('[INFO] TRANSFORMED THE DATA!!')
    print(f'data shape: {data_val.shape}')

    ## save pickle##
    l_value = []
    for l_drg, l_pat in itertools.product(drug_num, patient_num):
        l_ = data_val.query(f'(drug_id == {l_drg}) & (patient_id == {l_pat})').groupby(
            'day')['value'].median().reset_index(drop=True)
        l_value.append(l_)
        #b = pd.DataFrame(l_value, index = pat_list)
        #l_df.append(b)
    l_value_ = list(np.array_split(l_value, 7))  # reshape the list
    l_gr = []  # list of tumor growth rates
    for i in drug_num:
        gr_ = pd.DataFrame(l_value_[i-1], index=patient_num, columns=day_list)
        l_gr.append(gr_)
    with open(f'{savepath}.pickle', mode='wb') as f:
        pickle.dump(l_gr, f)
        f.close
    ## save csv ##
    data_val.to_csv(f'{savepath}.csv')
    return data_val, patient_num, drug_num


def cal_pvalue(data, patient_num, drug_num):
    # make list
    drug_name_list = data["drug_name"].drop_duplicates().tolist()
    day_list = data["day"].drop_duplicates().tolist()
    day_list.remove(0)  # drop '0' day

    # prepare dataset to calculate p-value
    list_t = []
    for pat, day, drg in itertools.product(patient_num, day_list, drug_num):
        v = data[(data['patient_id'] == pat) & (data['day'] == day)
                 & (data['drug_id'] == drg)]['value']
        list_t.append(v)
        a = np.array(list_t)
    data_reshape = a.reshape([80, 7, 6])

    # calculate p-value
    list_p = []
    for m in range(0, 80):
        for l in range(1, 7):  # t-test
            p = stats.ttest_rel(data_reshape[m][0], data_reshape[m][l])[1]
            list_p.append(p)
    data_pval_p = pd.DataFrame(np.reshape(
        list_p, (80, 6)), columns=drug_name_list[1:7])

    # make dataframe drug_name
    data_pval_f = data[["patient_id", "day"]].drop_duplicates()
    data_pval_f = data_pval_f[data_pval_f['day'] != 0].reset_index(drop=True)
    data_pval = pd.concat([data_pval_f, data_pval_p], axis=1)
    return data_pval


def fig_effects_per_patient(data, patient_num, estimator, ci):
    # option
    '''
    ## Args ##
    data : transformed InVivo data
    patient_num : patient number [1~10]
    estimator : (optional) Statistical function to estimate [mean or median]
    ci : (optional) Size of confidence intervals to draw around estimated values [float or “sd”]
    '''
    ## Args [estimator] ##
    if estimator == 'median':
        e = np.median
    else:
        estimator = 'mean'
        e = np.mean

    ## Args [ci] ##
    if type(ci) is int:
        ci = ci
    else:
        ci = 'sd'

    ## Args [data, patient_num] ##
    if type(patient_num) is int:
        data = data[data["patient_id"] == patient_num]
        patient_num = patient_num
        savepath = f'resultA_TumorGrowthRate/plot_{estimator}/fig_{estimator}_ci_{ci}_patient_{patient_num}.png'
    else:
        data = data
        patient_num = 'All [1~10]'
        savepath = f'resultA_TumorGrowthRate/plot_{estimator}/fig_{estimator}_ci_{ci}_patient_all.png'

    ## catplot ##
    # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
    colorlist = ["#999999", # gray
                "#9BC99B","#3FBF3F", # green
                "#9BA7C9", "#3F5FBF", # blue
                "#D88C8C", "#BF3F3F" # red
                ]
    sns.set_palette(colorlist, 7)

    sns.catplot(x='day', y='value', data=data, kind='point', hue='drug_name',
                estimator=e, ci=ci, height=8, aspect=2)
    # supplementary wire [pink]
    plt.hlines(1, -0.3, 8.3, '#D88C8C', linestyles='dashed')
    sns.set_style("darkgrid", {'grid.linestyle': '--'})
    plt.title(f'Patient {patient_num}, {estimator}, CI = {ci}', fontsize=20)
    ## save ##
    plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
    print(f'[SAVE]: {savepath}')


def fig_effects_per_drug(data, drug_num, estimator, ci):
    ## Args [estimator] ##
    if estimator == 'median':
        e = np.median
    else:
        estimator = 'mean'
        e = np.mean

    ## Args [ci] ##
    if type(ci) is int:
        ci = ci
    else:
        ci = 'sd'

    ## Args [data, drug_num] ##
    if type(drug_num) is int:
        data = data[data["drug_id"] == drug_num]
        drug_name = data[data["drug_id"] == drug_num]["drug_name"].iloc[1]
        savepath = f'resultA_TumorGrowthRate/plot_{estimator}/fig_{estimator}_ci_{ci}_drug_{drug_name}.png'
    else:
        data = data
        drug_name = 'All [1~10]'
        savepath = f'resultA_TumorGrowthRate/plot_{estimator}/fig_{estimator}_ci_{ci}_drug_all.png'

    sns.set_palette('Set3')
    sns.catplot(x='day', y='TumorGrowthRate', data=data, kind='point', hue='patient_id', estimator=e, ci=ci, height=8, aspect=2)
    #  Defaults are size=5, aspect=1
    plt.hlines(1, -0.3, 8.3, '#D88C8C', linestyles='dashed')  # 補助線
    sns.set_style("darkgrid", {'grid.linestyle': '--'})
    plt.title(f'Drug: {drug_name}, {estimator}, CI = {ci}', fontsize=20)

    plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
    print(f'[SAVE]: {savepath}')


def clustering(data):
    ## make list##
    drug_list = data["drug_id"].drop_duplicates().tolist()
    patient_list = data["patient_id"].drop_duplicates().tolist()
    day_list = data["day"].drop_duplicates().tolist()
    drug_index = data.drop_duplicates(
        'drug_id').iloc[:, 1:3].reset_index(drop=True)

    l_value = []
    for l_drug, l_patient in itertools.product(drug_list, patient_list):
        l_ = data.query(f'(drug_id == {l_drug}) & (patient_id == {l_patient})').groupby(
            'day')['value'].median().reset_index(drop=True)
        l_value.append(l_)
    l_value_ = list(np.array_split(l_value, 7))
    l_gr = []
    for i in drug_list:
        gr_ = pd.DataFrame(l_value_[i-1], index=patient_list, columns=day_list)
        l_gr.append(gr_)

    ## plot ##
    for i in range(len(drug_index)):
        # read data
        d = l_gr[i]
        labels = list(d.index)
        # clustering (method = ward)
        z = linkage(d, method='ward')
        # plot
        plt.figure(figsize=(30, 15))
        plt.rcParams["font.family"] = 'sans-serif'
        plt.rcParams['font.size'] = 30
        plt.title(
            f'Clustering Dendrogram : {drug_index["drug_name"][i]}', fontsize=50)
        plt.xlabel('CTOS_line', fontsize=40)
        plt.ylabel('Distance', fontsize=40)
        set_link_color_palette(['red', 'blue'])  # 2 cluster
        dendrogram(z, leaf_font_size=40, color_threshold=7.,
                   labels=labels, above_threshold_color='black')
        # save & show
        savepath = f'resultA_TumorGrowthRate/Clustering/fig_clustering_{drug_index["drug_name"][i]}.png'
        plt.savefig(savepath, dpi=300, format='png',
                    bbox_inches="tight")  # save
        print(f'[SAVE]: {savepath}')
    return


if __name__ == '__main__':
    ## dataload ##
    rawdata = dataload('data/InVivoData/InVivo_TumorGrowthRate.csv')
    ## transforming ##
    savepath = "data/data_TumorGrowthRate/data_TumorGrowthRate"
    TGR_data, patient_num, drug_num = data_transform(rawdata, savepath)
    print(f'[SAVE]: tumor growth rate data') # save format -> .csv, .pickle

    ## calculate p-value ##
    data_pval = cal_pvalue(TGR_data, patient_num, drug_num)
    data_pval.to_csv('data/data_TumorGrowthRate/pvalue_TumorGrowthRate.csv')
    print(f'[SAVE]: p-value data')

    ## Figure write [per patient] ##
    for l in patient_num:
        fig_effects_per_patient(TGR_data, l, 'mean', 'sd')  # mean plot
        fig_effects_per_patient(TGR_data, l, 'median', 'sd')  # median plot

    for d in drug_num:
        fig_effects_per_drug(TGR_data, d, 'mean', 'sd')
        fig_effects_per_drug(TGR_data, d, 'median', 'sd')

    # Figure write [all patient]
    fig_effects_per_patient(TGR_data, all, 'mean', 'sd')  # mean plot
    fig_effects_per_patient(TGR_data, all, 'median', 'sd')  # median plot

    ## clustering ##
    clustering(TGR_data)
