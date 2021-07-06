# Author: S. Inoue and yoshi
# Date: 9/18/2020
# Updated: 11/16/2020
# Project: CTOS
# Scropt: To generate dataset from raw microarray data for BN input
# Array dataset: CRC or SCNEC

import pandas as pd
import numpy as np
import glob
import os

# generate_data ---> input & sort & group
def generate_data(file):
    # read the data
    data = pd.read_table(str(file), sep='\t', header=9)
    print(f'input data matrix: {data.shape}')

    data_ = data.loc[:, ['ControlType',
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

    # sort
    data_ = data_[data_['ControlType'] == 0]
    data_ = data_[data_['gIsSaturated'] == 0]
    data_ = data_[data_['gIsFeatNonUnifOL'] == 0]
    data_ = data_[data_['gIsBGNonUnifOL'] == 0]
    data_ = data_[data_['gIsFeatPopnOL'] == 0]
    data_ = data_[data_['gIsBGPopnOL'] == 0]
    data_ = data_[data_['gIsWellAboveBG'] == 1]
    # data_ = data_[data_['gIsPosAndSignif'] == 0] ?
    print(f'sorted data matrix: {data_.shape}')

    # grouping
    generated_data_ = data_.groupby(
        'GeneName')['gProcessedSignal'].mean().sort_index()
    print(f'output data matrix: {generated_data_.shape}')

    return generated_data_


# extracting a path_name
def get_path_name(file):
    path_name_ = os.path.splitext(os.path.basename(file))[0]
    print(f'file name: {path_name_}')
    return path_name_


# concat allow duplication
def concat_outer(data):
    data_concat_outer_ = pd.concat(data, join='outer', axis=1)
    data_concat_outer_ = np.log2(data_concat_outer_ + 1)  # log2 transformation
    data_concat_outer_ = data_concat_outer_.sort_index()
    print(f'[OUTPUT] concat data matrix: {data_concat_outer_.shape}')
    return data_concat_outer_


# concat drop duplication
def concat_inner(data):
    data_concat_inner_ = pd.concat(data, join='inner', axis=1)
    data_concat_inner_ = np.log2(data_concat_inner_ + 1)  # log2 transformation
    data_concat_inner_ = data_concat_inner_.sort_index()
    print(f'[OUTPUT] final data matrix: {data_concat_inner_.shape}')
    return data_concat_inner_


# save file
def save_file(filename, data):
    data.to_csv(filename, mode='w', sep="\t")
    print(f"[SAVE] {filename}")


if __name__ == '__main__':
    file_list = sorted(glob.glob("../ArrayDataset/*.txt")) # set dataset path: CRC or SCNEC
    # generate data
    list = []
    for file in file_list:
        path_name = get_path_name(file)
        generated_data = generate_data(file).rename(path_name)
        list.append(generated_data)
    else:
        print('[INFO] Data generation completed')

    # concat
    # data_row = concat_outer(list)  # allow duplication
    data_final = concat_inner(list)  # drop duplication

    # save
    # save_file('data_row.txt', data_row)
    save_file('../data/dataset_final.txt', data_final)

