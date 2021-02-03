# Author: S.Inoue
# Date: 01/26/2021
# Updated: 01/26/2021
# Project: CTOS folfoliox project

# import module
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.model_selection import LeaveOneOut
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


# load the data
def data_load():
    # [LOAD] SelectedEdges Cetuximab_60mg (cxm60)
    file_1 = 'data/SelectedEdges_CorrCoef_rangeECv/SelectedEdges_CorrCoef0.6_rangeECv1.0_Cetuximab_60mg.txt'
    df_edge_cxm60 = pd.read_table(file_1, sep='\t', header=0)
    print(f'[LOAD]: {file_1}, input matrix: {df_edge_cxm60.shape}')

    # [LOAD] SelectedEdges Irinotecan_10mg (irn10)
    file_2 = 'data/SelectedEdges_CorrCoef_rangeECv/SelectedEdges_CorrCoef0.6_rangeECv1.0_Irinotecan_10mg.txt'
    df_edge_irn10 = pd.read_table(file_2, sep='\t', header=0)
    print(f'[LOAD]: {file_2}, input matrix: {df_edge_irn10.shape}')

    # [LOAD] SelectedEdges Oxaliplatin_10mg (oxa10)
    file_3 = 'data/SelectedEdges_CorrCoef_rangeECv/SelectedEdges_CorrCoef0.6_rangeECv1.0_Oxaliplatin_10mg.txt'
    df_edge_oxa10 = pd.read_table(file_3, sep='\t', header=0)
    print(f'[LOAD]: {file_3}, input matrix: {df_edge_oxa10.shape}')

    # [LOAD] original ECv data
    file_4 = 'data/ECv_network_result_0.15.txt'
    data_ecv_ = pd.read_table(file_4, sep='\t', header=0)
    print(f'[LOAD]: {file_4}, input matrix: {data_ecv_.shape}')

    return df_edge_cxm60, df_edge_irn10, df_edge_oxa10, data_ecv_

# reshape
def reshape(data_ecv_):
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

    print('[INFO] reshape the data')
    return data_ecv


# make input dataset
def preprocessing(data_ecv, df_edge_cxm60, df_edge_irn10, df_edge_oxa10):
    # drug list
    drug_list = [df_edge_cxm60, df_edge_irn10, df_edge_oxa10]
    # set the label
    label_cxm60 = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1]  # 2,10
    label_irn60 = [0, 0, 0, 1, 0, 1, 0, 0, 0, 1]  # 4,6,10
    label_oxa60 = [0, 0, 1, 1, 0, 0, 1, 0, 1, 1]  # 3,4,7,9,10
    label = [label_cxm60, label_irn60, label_oxa60]
    # 各薬剤で成形
    reshaped = []
    for i in range(len(drug_list)):
        # 検索対象となるペア列を作成
        data_ecv["pair"] = data_ecv["Parent"] + "-" + data_ecv["Child"]
        drug_list[i]["pair"] = drug_list[i]["Parent"] + \
            "-" + drug_list[i]["Child"]
        # data_ecvから必要な列の抽出
        selected_ecv = data_ecv[data_ecv['pair'].isin(list(drug_list[i]["pair"]))].drop(
            columns=['Parent', 'Child', 'pair']).T  # 不要な行を削除して転置
        # ラベル列を設置
        selected_ecv["label"] = label[i]
        reshaped.append(selected_ecv)
        # df_reshaped[0~2] [0]cxm60 [1]irn10 [2]oxa10

    return reshaped


# learning
def learning(df_reshaped):
    method = ['SVM', 'LDA']
    for m in range(len(method)):
        print(f'method : {method[m]}')
        for i in range(len(df_reshaped)):
            # drug list
            drug_name = ["Cetuximab_60mg", "Irinotecan_10mg", "Oxaliplatin_10mg"]
            print(f'drug : {drug_name[i]}')

            data = np.array(df_reshaped[i].drop(columns='label'))  # numpy行列に変換
            target = np.array(df_reshaped[i]["label"])  # numpy行列に変換

            loo = LeaveOneOut()  # LOOCVのインスタンス生成

            # インスタンスの選択
            if m == 0:
                clf = svm.SVC(gamma="scale")
            elif m == 1:
                clf = LinearDiscriminantAnalysis()

            entire_count = loo.get_n_splits(data)  # テスト回数取得
            correct_answer_count = 0  # 推定が正解だった数初期化

            print("[predict] [Answer]")
            # loo.split(data)で訓練データとテストデータを分割
            for train_index, test_index in loo.split(data):
                #print("TRAIN:", train_index, "TEST:", test_index)
                data_train, data_test = data[train_index], data[test_index]
                target_train, target_test = target[train_index], target[test_index]
                clf.fit(data_train, target_train)  # 学習させる
                result = clf.predict(data_test)  # テストデータからラベルを予測する
                print(f'    {result} {target_test}')  # 答えを表示
                if result == target_test:  # ラベルと元々のラベルが一致していれば+1
                    correct_answer_count += 1

            rate = (float(correct_answer_count) /
                    float(entire_count))  # 正解率を計算
            print(f'  Correct Answer Rate : {str(rate)}')  # 正解率を出力
    return


if __name__ == '__main__':
    # load the data
    df_edge_cxm60, df_edge_irn10, df_edge_oxa10, data_ecv_ = data_load()
    # reshape the data
    data_ecv = reshape(data_ecv_)
    # make input data
    df_reshaped = preprocessing(data_ecv, df_edge_cxm60, df_edge_irn10, df_edge_oxa10)
    # learning
    learning(df_reshaped)
