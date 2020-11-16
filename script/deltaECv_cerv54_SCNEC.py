# Author: yoshi
# Date: 11/05/2020
# Updated:
# Project: CTOS SCNEC
# Script: to calculate delta ECv

import pandas as pd
import matplotlib.pyplot as plt

def main():

    ecv = pd.read_table('../NetworkFiles/SCNEC/ECv_network_result_0.10_processing.txt', sep='\t', header=0)
    # SCNEC samples: cerv5,8,9,21,39,46,51,54,59,60,61,TY-YIK

    ## deltaECv = 54 vs XX
    ecv['deltaECv_54vs5'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv5:5'])
    ecv['deltaECv_54vs8'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv8:1'])
    ecv['deltaECv_54vs9'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv9:6'])
    ecv['deltaECv_54vs21'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv21:2'])
    ecv['deltaECv_54vs39'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv39:3'])
    ecv['deltaECv_54vs46'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv46:4'])
    ecv['deltaECv_54vs51'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv51:7'])
    ecv['deltaECv_54vs59'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv59:10'])
    ecv['deltaECv_54vs60'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv60:11'])
    ecv['deltaECv_54vs61'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:cerv61:12'])
    ecv['deltaECv_54vsTYYIK'] = abs(ecv['ECv:cerv54:9'] - ecv['ECv:TY-YIK:8'])


    ## plot deltaECv ditribution
    delta_54vs5 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs5']]
    plt.hist(delta_54vs5['deltaECv_54vs5'], log=True, color='rosybrown', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename1 = '../data/SCNEC_deltaECv_histogram_54vs5.png'
    print(f'[SAVE]: {filename1}')
    plt.savefig(filename1, dpi=300, format='png')
    plt.clf()

    delta_54vs8 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs8']]
    plt.hist(delta_54vs8['deltaECv_54vs8'], log=True, color='darksalmon', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename2 = '../data/SCNEC_deltaECv_histogram_54vs8.png'
    print(f'[SAVE]: {filename2}')
    plt.savefig(filename2, dpi=300, format='png')
    plt.clf()

    delta_54vs9 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs9']]
    plt.hist(delta_54vs9['deltaECv_54vs9'], log=True, color='sandybrown', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename3 = '../data/SCNEC_deltaECv_histogram_54vs9.png'
    print(f'[SAVE]: {filename3}')
    plt.savefig(filename3, dpi=300, format='png')
    plt.clf()

    delta_54vs21 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs21']]
    plt.hist(delta_54vs21['deltaECv_54vs21'], log=True, color='darkkhaki', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename4 = '../data/SCNEC_deltaECv_histogram_54vs21.png'
    print(f'[SAVE]: {filename4}')
    plt.savefig(filename4, dpi=300, format='png')
    plt.clf()

    delta_54vs39 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs39']]
    plt.hist(delta_54vs39['deltaECv_54vs39'], log=True, color='olivedrab', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename5 = '../data/SCNEC_deltaECv_histogram_54vs39.png'
    print(f'[SAVE]: {filename5}')
    plt.savefig(filename5, dpi=300, format='png')
    plt.clf()

    delta_54vs46 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs46']]
    plt.hist(delta_54vs46['deltaECv_54vs46'], log=True, color='chartreuse', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename6 = '../data/SCNEC_deltaECv_histogram_54vs46.png'
    print(f'[SAVE]: {filename6}')
    plt.savefig(filename6, dpi=300, format='png')
    plt.clf()

    delta_54vs51 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs51']]
    plt.hist(delta_54vs51['deltaECv_54vs51'], log=True, color='seagreen', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename7 = '../data/SCNEC_deltaECv_histogram_54vs51.png'
    print(f'[SAVE]: {filename7}')
    plt.savefig(filename7, dpi=300, format='png')
    plt.clf()

    delta_54vs59 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs59']]
    plt.hist(delta_54vs59['deltaECv_54vs59'], log=True, color='darkcyan', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename8 = '../data/SCNEC_deltaECv_histogram_54vs59.png'
    print(f'[SAVE]: {filename8}')
    plt.savefig(filename8, dpi=300, format='png')
    plt.clf()

    delta_54vs60 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs60']]
    plt.hist(delta_54vs60['deltaECv_54vs60'], log=True, color='slategray', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename9 = '../data/SCNEC_deltaECv_histogram_54vs60.png'
    print(f'[SAVE]: {filename9}')
    plt.savefig(filename9, dpi=300, format='png')
    plt.clf()

    delta_54vs61 = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vs61']]
    plt.hist(delta_54vs61['deltaECv_54vs61'], log=True, color='darkturquoise', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename10 = '../data/SCNEC_deltaECv_histogram_54vs61.png'
    print(f'[SAVE]: {filename10}')
    plt.savefig(filename10, dpi=300, format='png')
    plt.clf()

    delta_54vsTYYIK = ecv.loc[:,['Parent', 'Child', 'deltaECv_54vsTYYIK']]
    plt.hist(delta_54vsTYYIK['deltaECv_54vsTYYIK'], log=True, color='plum', alpha=0.6, edgecolor='k', lw=0.8, bins=40)
    plt.title('ΔECv histogram')
    plt.xlabel('ΔECv')
    plt.ylabel('frequency')
    filename11 = '../data/SCNEC_deltaECv_histogram_54vsTYYIK.png'
    print(f'[SAVE]: {filename11}')
    plt.savefig(filename11, dpi=300, format='png')
    plt.clf()


    ## deltaECv threshold
    deltaECv_th = 1

    delta_54vs5_th = delta_54vs5[delta_54vs5['deltaECv_54vs5'] >= deltaECv_th]
    filename101 = '../data/SCNEC_deltaECv_table_54vs5_th1.txt'
    print(f'[SAVE]: {filename101}')
    with open(filename101, 'w') as f:
        delta_54vs5_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs8_th = delta_54vs5[delta_54vs8['deltaECv_54vs8'] >= deltaECv_th]
    filename102 = '../data/SCNEC_deltaECv_table_54vs8_th1.txt'
    print(f'[SAVE]: {filename102}')
    with open(filename102, 'w') as f:
        delta_54vs8_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs9_th = delta_54vs9[delta_54vs9['deltaECv_54vs9'] >= deltaECv_th]
    filename103 = '../data/SCNEC_deltaECv_table_54vs9_th1.txt'
    print(f'[SAVE]: {filename103}')
    with open(filename103, 'w') as f:
        delta_54vs9_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs21_th = delta_54vs21[delta_54vs21['deltaECv_54vs21'] >= deltaECv_th]
    filename104 = '../data/SCNEC_deltaECv_table_54vs21_th1.txt'
    print(f'[SAVE]: {filename104}')
    with open(filename104, 'w') as f:
        delta_54vs21_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs39_th = delta_54vs39[delta_54vs39['deltaECv_54vs39'] >= deltaECv_th]
    filename105 = '../data/SCNEC_deltaECv_table_54vs39_th1.txt'
    print(f'[SAVE]: {filename105}')
    with open(filename105, 'w') as f:
        delta_54vs39_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs46_th = delta_54vs46[delta_54vs46['deltaECv_54vs46'] >= deltaECv_th]
    filename106 = '../data/SCNEC_deltaECv_table_54vs46_th1.txt'
    print(f'[SAVE]: {filename106}')
    with open(filename106, 'w') as f:
        delta_54vs46_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs51_th = delta_54vs51[delta_54vs51['deltaECv_54vs51'] >= deltaECv_th]
    filename107 = '../data/SCNEC_deltaECv_table_54vs51_th1.txt'
    print(f'[SAVE]: {filename107}')
    with open(filename107, 'w') as f:
        delta_54vs51_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs59_th = delta_54vs59[delta_54vs59['deltaECv_54vs59'] >= deltaECv_th]
    filename108 = '../data/SCNEC_deltaECv_table_54vs59_th1.txt'
    print(f'[SAVE]: {filename108}')
    with open(filename108, 'w') as f:
        delta_54vs59_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs60_th = delta_54vs60[delta_54vs60['deltaECv_54vs60'] >= deltaECv_th]
    filename109 = '../data/SCNEC_deltaECv_table_54vs60_th1.txt'
    print(f'[SAVE]: {filename109}')
    with open(filename109, 'w') as f:
        delta_54vs60_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vs61_th = delta_54vs61[delta_54vs61['deltaECv_54vs61'] >= deltaECv_th]
    filename110 = '../data/SCNEC_deltaECv_table_54vs61_th1.txt'
    print(f'[SAVE]: {filename110}')
    with open(filename110, 'w') as f:
        delta_54vs61_th.to_csv(f, sep='\t', header=True, index=False)

    delta_54vsTYYIK_th = delta_54vsTYYIK[delta_54vsTYYIK['deltaECv_54vsTYYIK'] >= deltaECv_th]
    filename111 = '../data/SCNEC_deltaECv_table_54vsTYYIK_th1.txt'
    print(f'[SAVE]: {filename111}')
    with open(filename111, 'w') as f:
        delta_54vsTYYIK_th.to_csv(f, sep='\t', header=True, index=False)





if __name__ == '__main__':
    main()

