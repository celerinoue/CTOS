# Author: S.Inoue
# Date: 02/01/2021
# Updated: 02/01/2021
# Project: CTOS folfoliox project


#%%
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

# %%
def learning(df_all):
    # #############################################################################
    # name list
    name = []
    edge_name = ["CorrCoef0.6_rangeECv1.0", "CorrCoef0.6_rangeECv1.0_forGeneExp",
                 "CorrCoef0.7_rangeGeneExp1.0", "CorrCoef0.8_rangeGeneExp1.0"]
    drug_name = ["Cetuximab_60mg", "Irinotecan_10mg", "Oxaliplatin_10mg"]
    for edge, drug in itertools.product(edge_name, drug_name):
        name.append([edge, drug])


    for i in range(len(df_all)):
        # #############################################################################
        # Data Name
        print(f'edge: {name[i][0]}, drug: {name[i][1]}')

        # #############################################################################
        # Data IO and generation
        X = np.array(df_all[i].drop(columns='label'))  # Explanatory variables
        y = np.array(df_all[i]["label"])  # Dependent variable
        print(f'[INFO] data shape: X {X.shape}, y {y.shape}')
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
            all_y.append(y[test]) # true label
            all_probs.append(clf.predict_proba(X[test])[:, 1]) # probability
            all_preds.append(clf.predict(X[test]))  # predict value

        all_y = np.ravel(all_y)
        all_probs = np.ravel(all_probs)
        all_preds = np.ravel(all_preds)

        # count correct_num
        correct_ = []
        for n in range(len(all_y)):
            correct_.append(all_preds[n] == all_y[n])
            correct_rate = correct_.count(True) / 10

        # plot
        fpr, tpr, thresholds = roc_curve(all_y, all_probs)
        roc_auc = auc(fpr, tpr)

        print(f"TRUE: {all_y}, PRED: {all_preds}, prob: {all_probs}, correct rate: {correct_rate}, AUC:{roc_auc}", fpr, tpr)


        plt.figure(1, figsize=(12, 6))
        plt.plot(fpr, tpr, lw=2, alpha=0.5,
                label='LOOCV ROC (AUC = %0.2f)' % (roc_auc))
        plt.plot([0, 1], [0, 1], linestyle='--', lw=2,
                color='k', label='Chance level', alpha=.2)

        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(
            f'CTOS classification ROC Curves [method:Linear Discriminant Analysis][Edge:{name[i][0]}][Drug: {name[i][1]}]')
        plt.legend(loc="lower right")
        plt.grid()

        #savepath = f'result/fig/IDEA1_4/1_classification_ROC_curves/fig_ROC_{name[i][0]}_{name[i][1]}.png'
        #plt.savefig(savepath, dpi=300, format='png', bbox_inches="tight")
        #print(f'[SAVE]: {savepath}')


    return

learning(df_all)




if __name__ == '__main__':
    # load the data
    df_all = data_load()
    # predict & plot ROC curves
    learning(df_all)
