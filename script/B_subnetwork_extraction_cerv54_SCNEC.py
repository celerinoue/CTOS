# Author: yoshi
# Date: 11/10/2020
# Updated:
# Project: CTOS SCNEC
# Script: To extract core subnetwork for cerv54

import pandas as pd

def main():

    vs5=pd.read_table('SCNEC_deltaECv_table_54vs5_th1.txt',sep='\t',header=0)
    vs8=pd.read_table('SCNEC_deltaECv_table_54vs8_th1.txt',sep='\t',header=0)
    vs9=pd.read_table('SCNEC_deltaECv_table_54vs9_th1.txt',sep='\t',header=0)
    vs21=pd.read_table('SCNEC_deltaECv_table_54vs21_th1.txt',sep='\t',header=0)
    vs39=pd.read_table('SCNEC_deltaECv_table_54vs39_th1.txt',sep='\t',header=0)
    vs46=pd.read_table('SCNEC_deltaECv_table_54vs46_th1.txt',sep='\t',header=0)
    vs51=pd.read_table('SCNEC_deltaECv_table_54vs51_th1.txt',sep='\t',header=0)
    vs59=pd.read_table('SCNEC_deltaECv_table_54vs59_th1.txt',sep='\t',header=0)
    vs60=pd.read_table('SCNEC_deltaECv_table_54vs60_th1.txt',sep='\t',header=0)
    vs61=pd.read_table('SCNEC_deltaECv_table_54vs61_th1.txt',sep='\t',header=0)

    vs5_list=[]
    vs8_list=[]
    vs9_list=[]
    vs21_list=[]
    vs39_list=[]
    vs46_list=[]
    vs51_list=[]
    vs59_list=[]
    vs60_list=[]
    vs61_list=[]

    for k,l in zip(vs5['Parent'],vs5['Child']):
        pair=k+'_'+l
        vs5_list.append(pair)

    for k,l in zip(vs8['Parent'],vs8['Child']):
        pair=k+'_'+l
        vs8_list.append(pair)

    for k,l in zip(vs9['Parent'],vs9['Child']):
        pair=k+'_'+l
        vs9_list.append(pair)

    for k,l in zip(vs21['Parent'],vs21['Child']):
        pair=k+'_'+l
        vs21_list.append(pair)

    for k,l in zip(vs39['Parent'],vs39['Child']):
        pair=k+'_'+l
        vs39_list.append(pair)

    for k,l in zip(vs46['Parent'],vs46['Child']):
        pair=k+'_'+l
        vs46_list.append(pair)

    for k,l in zip(vs51['Parent'],vs51['Child']):
        pair=k+'_'+l
        vs51_list.append(pair)

    for k,l in zip(vs59['Parent'],vs59['Child']):
        pair=k+'_'+l
        vs59_list.append(pair)

    for k,l in zip(vs60['Parent'],vs60['Child']):
        pair=k+'_'+l
        vs60_list.append(pair)

    for k,l in zip(vs61['Parent'],vs61['Child']):
        pair=k+'_'+l
        vs61_list.append(pair)


    shared_edges = (set(vs5_list) &
                    set(vs8_list) &
                    set(vs9_list) &
                    set(vs21_list) &
                    set(vs39_list) &
                    set(vs46_list) &
                    set(vs59_list) &
                    set(vs60_list) &
                    set(vs61_list))

    print(shared_edges)


if __name__ == '__main__':
    main()
