#%%
import pandas as pd
import numpy as np
import glob
import os
import GEOparse


#%%
# load
def load_data(file):
    ## load data ##
    data_raw = pd.read_table(str(file), sep='\t')
    print(f'# input data matrix: {data_raw.shape}')

    return data_raw


if __name__ == '__main__':
    data_raw = load_data(
        '/Users/celerinoue/0_res/CTOS/data_FOLFOXRF/GSE28702_family.soft')
    data_raw
# %%

gse = GEOparse.get_GEO(filepath="data_FOLFOXRF/GSE28702_family.soft")
platform = gse.gpls['GPL570'].columns

# %%
