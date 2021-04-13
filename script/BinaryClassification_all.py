# Author: S.Inoue
# Date: 02/01/2021
# Updated: 03/18/2021
# Project: CTOS folfoliox project


# import module
import numpy as np
import pandas as pd
import glob
import os
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import LeaveOneOut
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_curve, auc


# load selected CTOS set matrix
def data_load():
    df_all = []
    for i in sorted(glob.glob("data_BinaryClassification/all/*.txt")):
        print(f'[LOAD] {os.path.basename(i)}')
        df = pd.read_table(i, sep='\t', index_col=0)
        df_all.append(df)
    return df_all


def learning(df_all):
    # name list
    name = []
    edge_name = ["CorrCoef0.6_rangeECv1.0", "CorrCoef0.6_rangeECv1.0_forGeneExp",
                 "CorrCoef0.7_rangeGeneExp1.0", "CorrCoef0.8_rangeGeneExp1.0"]
    drug_name = ["Cetuximab_60mg", "Irinotecan_10mg", "Oxaliplatin_10mg"]
    for edge, drug in itertools.product(edge_name, drug_name):
        name.append([edge, drug])

    # learning
    print('[INFO] binary classification ===============')
    all_auc, all_accuracy = [], []

    for i in range(len(df_all)):
        print(f'[INFO] data_num {i}')
        # Data Name
        print(f'# edge: {name[i][0]}, drug: {name[i][1]}')
        # Data IO and generation
        X = np.array(df_all[i].drop(columns='label'))  # Explanatory variables
        y = np.array(df_all[i]["label"])  # Dependent variable
        #print(f'# data shape: X {X.shape}, y {y.shape}')

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
            all_y.append(y[test]) # true label
            all_probs.append(clf.predict_proba(X[test])[:, 1]) # probability
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

    return all_accuracy

def plot(all_accuracy):
    # y
    y = [np.array(all_accuracy[0:3]),  # acc_th0.6_ECv
         np.array(all_accuracy[3:6]),  # acc_th0.6_GeneExp
         np.array(all_accuracy[6:9]),  # acc_th0.7_GeneExp
         np.array(all_accuracy[9:12]),  # acc_th0.8_GeneExp
         ]
    # x
    x = list(range(1, 4))  # x range

    # title, legend
    node = ['ECv', 'GeneExpression']
    num_node = [0, 1, 1, 1]
    drug = ['Cetuximab_60mg', 'Irinotecan_10mg', 'Oxaliplatin_10mg']
    th = ['0.6', '0.6', '0.7', '0.8']

    # color
    # (Cetuximab, Oxaliplatin, Irinotecan) = (green, blue, red)
    colorlist = ["#9BC99B", "#9BA7C9", "#D88C8C"]

    for i in range(len(y)):
        sns.set()
        sns.set_style("white")

        fig = plt.figure(figsize=(12, 9))
        ax = sns.barplot(x=x, y=y[i], palette=colorlist)
        ax.set_ylim([0, 1])
        ax.set_title(
            f'Accuracy of Binary Classification [Data = {node[num_node[i]]}, CorrCoef(th = {th[i]})]')
        ax.set_xlabel(f'drug name')
        ax.set_xticklabels(drug)
        ax.set_ylabel('Accuracy')

        savepath = f'resultB_BinaryClassification/Classification_result/fig_Classification_result_{node[num_node[i]]}_th{th[i]}_.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')

    return


if __name__ == '__main__':
    # load the data
    df_all = data_load()
    # predict
    all_accuracy = learning(df_all)
    # plot
    plot(all_accuracy)
