# Author: S.Inoue
# Date: 02/01/2021
# Updated: 02/01/2021
# Project: CTOS folfoliox project

# import module
import numpy as np
import pandas as pd
import glob
import os
import itertools
import matplotlib.pyplot as plt
from sklearn.model_selection import LeaveOneOut
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import roc_curve, auc


# load selected CTOS set matrix
def data_load():
    df_all = []
    for i in sorted(glob.glob("result/txt/IDEA1_4/SelectedCTOSset/*.txt")):
        print(f'[LOAD] {os.path.basename(i)}')
        df = pd.read_table(i, sep='\t', index_col=0)
        df_all.append(df)
    return df_all


def learning(df_all):
    # #############################################################################
    # name list
    name = []
    edge_name = ["CorrCoef0.6_rangeECv1.0", "CorrCoef0.6_rangeECv1.0_forGeneExp",
                 "CorrCoef0.7_rangeGeneExp1.0", "CorrCoef0.8_rangeGeneExp1.0"]
    drug_name = ["Cetuximab_60mg", "Irinotecan_10mg", "Oxaliplatin_10mg"]
    for e, d in itertools.product(edge_name, drug_name):
        name.append([e, d])


    for i in range(len(df_all)):
        # #############################################################################
        # Data Name
        print(f'edge: {name[i][0]}, drug: {name[i][1]}')

        # #############################################################################
        # Data IO and generation
        X = np.array(df_all[i].drop(columns='label'))  # Explanatory variables
        y = np.array(df_all[i]["label"])  # Dependent variable

        # #############################################################################
        # Classification and ROC analysis

        clf = LinearDiscriminantAnalysis(solver='svd',
                                         #shrinkage='auto',
                                         tol=1.0e-4)
        cv = LeaveOneOut()

        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        # Run classifier with LeaveOneOut and plot ROC curves
        fig = plt.figure(figsize=[12, 12])
        ax = fig.add_subplot(111, aspect='equal')
        print(name[i])

        for train, test in cv.split(X):
            #print("TRAIN:", train, "TEST:", test)
            clf.fit(X[train], y[train])
            pred = clf.predict(X)  # predict values
            true = y  # true values

            # count correct_num
            correct_ = []
            for n in range(len(pred)):
                correct_.append(pred[n] == y[n])
            correct_rate = correct_.count(True) / 10
            #print(f'TRUE: {true}, PRED: {pred}, correct rate: {correct_rate}')

            fpr, tpr, t = roc_curve(true, pred)
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
            plt.plot(fpr, tpr, lw=2, alpha=0.2,
                     label='ROC Leave-One-Out %d (AUC = %0.2f)' % (i, roc_auc))

        plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='black')
        mean_tpr = np.mean(tprs, axis=0)
        mean_auc = auc(mean_fpr, mean_tpr)
        plt.plot(mean_fpr, mean_tpr, color='blue',
                label=r'Mean ROC (AUC = %0.2f )' % (mean_auc), lw=2, alpha=1)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'CTOS classification ROC Curves [method:Linear Discriminant Analysis][Edge:{name[i][0]}][Drug: {name[i][1]}]')
        plt.legend(loc="lower right")

        savepath = f'result/fig/IDEA1_4/1_classification_ROC_curves/fig_ROC_{name[i][0]}_{name[i][1]}.png'
        plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        print(f'[SAVE]: {savepath}')
    return


if __name__ == '__main__':
    # load the data
    df_all = data_load()
    # predict & plot ROC curves
    learning(df_all)
