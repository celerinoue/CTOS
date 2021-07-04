# Author: S. Inoue and yoshi
# Date: 9/18/2020
# Updated: 03/18/2021
# Project: CTOS
# Scropt: To generate dataset from raw microarray data for BN input
# Array dataset: CTOS CRC

import pandas as pd
import numpy as np
import glob
import os

# feature extract ---> input & sort & grouping
def feature_extract(file):
    ## load data ##
    data_ = pd.read_table(str(file), sep='\t', header=9)
    print(f'# input data matrix: {data_.shape}')

    data = data_.loc[:, ['ControlType',
                         'ProbeName',
                         'GeneName',
                         'gProcessedSignal',
                         'gIsWellAboveBG',
                         'gIsSaturated',
                         'gIsFeatNonUnifOL',
                         'gIsBGNonUnifOL',
                         'gIsFeatPopnOL',
                         'gIsBGPopnOL',
                         'gIsPosAndSignif']]

    ## sort ##
    data = data[data['ControlType'] == 0]
    data = data[data['gIsSaturated'] == 0]
    data = data[data['gIsFeatNonUnifOL'] == 0]
    data = data[data['gIsBGNonUnifOL'] == 0]
    data = data[data['gIsFeatPopnOL'] == 0]
    data = data[data['gIsBGPopnOL'] == 0]
    data = data[data['gIsWellAboveBG'] == 1]
    # data = data[data['gIsPosAndSignif'] == 0] ?

    ## grouping ##
    extracted_data = data.groupby('GeneName')['gProcessedSignal'].mean().sort_index()
    print(f'# sorted data matrix: {extracted_data.shape}')
    return extracted_data

# concat allow duplication
def concat_outer(list_data):
    data_concat_outer_ = pd.concat(list_data, join='outer', axis=1)
    data_concat_outer_ = np.log2(data_concat_outer_ + 1)  # log2 transformation
    data_concat_outer = data_concat_outer_.sort_index()
    print(f'# extracted data shape: {data_concat_outer.shape}')
    return data_concat_outer


# concat drop duplication
def concat_inner(list_data):
    data_concat_inner_ = pd.concat(list_data, join='inner', axis=1)
    data_concat_inner_ = np.log2(data_concat_inner_ + 1)  # log2 transformation
    data_concat_inner = data_concat_inner_.sort_index()
    print(f'# extracted data shape: {data_concat_inner.shape}')
    return data_concat_inner


# save file
def save_file(filename, data):
    data.to_csv(filename, mode='w', sep="\t")
    print(f"[SAVE] {filename}")


if __name__ == '__main__':
    ## load CRC data ##
    file_list = sorted(glob.glob('data/ArrayDataset/CRC/*.txt'))
    ## processing data ##
    list_extracted_data = []
    for file in file_list:
        path_name = os.path.splitext(os.path.basename(file))[0]
        print(f'# file name: {path_name}')
        f = feature_extract(file).rename(path_name)
        list_extracted_data.append(f)
    else:
        print('[INFO] feature extraction completed')

    ## concat ##
    # FeatureExtractedMatrix_raw = concat_outer(list_extracted_data)  # allow gene duplication
    FeatureExtractedMatrix = concat_inner(list_extracted_data)  # drop gene duplication

    ## save ##
    # save_file('BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset_raw.txt', FeatureExtractedMatrix_raw)
    save_file('BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt', FeatureExtractedMatrix)
