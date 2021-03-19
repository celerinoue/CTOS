# Author: S.Inoue
# Date: 2/22/2021
# Updated: 2/22/2021
# Project: CTOS folfoli folfox
# Script: Exract connected component from the network file generated through INGOR


# import module
import numpy as np
import pandas as pd
import glob
import os
import itertools
from natsort import natsorted
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import LeaveOneOut
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_curve, auc

os.chdir('/Users/celerinoue/0_res/CTOS/')

# load selected CTOS set matrix
def data_load():
    df_all, list_df_name = [], []
    for i in natsorted(glob.glob("data_BinaryClassification/SelectedComponent/*.txt")):
        print(f'[LOAD] {os.path.basename(i)}')
        df = pd.read_table(i, sep='\t', index_col=0)
        df_all.append(df)
        list_df_name.append([os.path.basename(i).split('_')[1], # 'GeneExp' or 'ECv'
                             os.path.basename(i).split('_')[2], # drug
                             os.path.basename(i).split('_')[5].split('.')[0]  # rank
                             ])
    return df_all, list_df_name


def learning(df_all, list_df_name):
    print('[INFO] binary classification ===============')

    all_auc, all_accuracy = [], []
    for i in range(len(df_all)):
        # Data Name
        print(f'[INFO] data_num {i}')
        print(f'# type: {list_df_name[i][0]}, drug: {list_df_name[i][1]}, rank: {list_df_name[i][2]}')

        # #############################################################################
        # Data IO and generation
        X = np.array(df_all[i].drop(columns='label'))  # Explanatory variables
        y = np.array(df_all[i]["label"])  # Dependent variable
        #print(f'# data shape: X {X.shape}, y {y.shape}')

        # #############################################################################
        # Classification and ROC analysis

        clf = LinearDiscriminantAnalysis(solver='svd',
                                         #shrinkage='auto',
                                         tol=1.0e-4)
        cv = LeaveOneOut()

        all_y = []
        all_probs = []
        all_preds = []

        for train, test in cv.split(X, y):
            clf.fit(X[train], y[train])
            all_y.append(y[test])  # true label
            all_probs.append(clf.predict_proba(X[test])[:, 1])  # probability
            all_preds.append(clf.predict(X[test]))  # predict value

        all_y = np.ravel(all_y)
        all_probs = np.ravel(all_probs)
        all_preds = np.ravel(all_preds)

        # count accuracy, auc
        correct_ = []
        for n in range(len(all_y)):
            correct_.append(all_preds[n] == all_y[n])
            accuracy = correct_.count(True) / 10

        fpr, tpr, thresholds = roc_curve(all_y, all_probs)
        roc_auc = auc(fpr, tpr)
        all_auc.append(roc_auc)
        all_accuracy.append(accuracy)

        #print(f"# TRUE: {all_y}, PRED: {all_preds}, prob: {all_probs}")
        print(f"# accuracy: {accuracy}, AUC:{roc_auc}")

    return all_accuracy, all_auc


def plot(all_accuracy, all_auc):
    # y
    y = [np.array(all_accuracy[0:14]),  # acc_ecv_cxm
        np.array(all_accuracy[14:20]),  # acc_ecv_irn
        np.array(all_accuracy[20:30]),  # acc_ecv_oxa
        np.array(all_accuracy[30:44]),  # acc_exp_cxm
        np.array(all_accuracy[44:50]),  # acc_exp_irn
        np.array(all_accuracy[50:60]),  # acc_exp_oxa
        np.array(all_auc[0:14]),  # auc_ecv_cxm
        np.array(all_auc[14:20]),  # auc_ecv_irn
        np.array(all_auc[20:30]),  # auc_ecv_oxa
        np.array(all_auc[30:44]),  # auc_exp_cxm
        np.array(all_auc[44:50]),  # auc_exp_irn
        np.array(all_auc[50:60])  # auc_exp_oxa
        ]
    # x
    x = [list(range(1, 15)), # cxm
        list(range(1, 7)), # irn
        list(range(1, 11))  # oxa
        ]
    # title
    edge = ['ECv', 'GeneExpression']
    num_edge = [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]
    method = ['Accuracy', 'AUC']
    num_method = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
    drug = ['Cetuximab_60mg', 'Irinotecan_10mg', 'Oxaliplatin_10mg']

    # color
    # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
    #colorlist = ["YlGn_r", "YlGnBu_r", "YlOrRd_r"]
    colorlist = ["Greens", "Blues", "Reds"]
    num_col = [1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3]

    for i in range(len(y)):
        sns.set()
        sns.set_style("white")
        sns.set_palette(colorlist[i % 3], num_col[i])

        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(1, 1, 1)
        ax.bar(x[i % 3], y[i], tick_label=x[i % 3])
        ax.set_ylim([0, 1])
        ax.set_title(f'{method[num_method[i]]} of binary classification [drug = {drug[i%3]}, edge = {edge[num_edge[i]]}, component length >= 4]')
        ax.set_xlabel(f'Ranking of Component Size [th = 0.6]')
        ax.set_ylabel(f'{method[num_method[i]]}')

        savepath = f'resultC_BinaryClassification_SelectedComponent/Classification_result/fig_Classification_result_{method[num_method[i]]}_{edge[num_edge[i]]}_{drug[i%3]}_.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')

    return


if __name__ == '__main__':
    # load data
    df_all, list_df_name = data_load()
    # classification
    all_accuracy, all_auc = learning(df_all, list_df_name)
    # save fig
    plot(all_accuracy, all_auc)
