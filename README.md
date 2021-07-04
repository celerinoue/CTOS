# CTOS_CRC

## A. 腫瘍増殖率のチェック, ラベリング
```
script/A_data_processing_TumorGrowthRate.py
```
1. データ前処理
    - `RawData/InVivoData/腫瘍径推移まとめ.xls`　から `RawData/InVivoData/InVivo_TumorGrowthRate.csv` を手動で抽出.
    - [output] 腫瘍径推移まとめ `data/data_TumorGrowthRate/data_TumorGrowthRate.csv`
    - [output] 腫瘍径推移 p-value 変化に有意な差あるかを確認 `data/data_TumorGrowthRate/pvalue_TumorGrowthRate.csv`

2. 可視化
    - [fig] 薬剤ごとの腫瘍径推移plot (マウスサンプルの平均値) `resultA_TumorGrowthRate/plot_mean/` 
    - [fig] 薬剤ごとの腫瘍径推移plot (マウスサンプルの中央値) `resultA_TumorGrowthRate/plot_median/` 
    - [fig] 腫瘍径推移 クラスタリング結果`resultA_TumorGrowthRate/Clustering/`
        - クラスタリング結果をもとにラベリング


## B. Bayesian Network推定

1. データ前処理
```
script/B_Preprocessing_BNinput_CTOS.py
```
- BN input matrix `BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt`


2. 可視化
```
script/B_data_processing.py
```
- 元データの分布を確認
    - [fig] ヒストグラム
    - [fig] PCA
    - [fig] violin plot

```
script/B_Visualization_CorrCoef_ECv.py
```
- ECvの分布を確認
    - [input] BN result (=ECvデータ) `BayesianNetworkEstimation/CRC/ECv_network_result_0.15.txt`
    - [output] ECvエッジと腫瘍増殖率の相関係数のリスト `data/data_CorrCoef/ECv/`
        - `Parent	Child	CorrelationCoefficients	range_ecv`
    - [fig] ECvのrangeの分布 `resultB_CorrCoef/range_ECv/dist_CorrCoef_ECv/`
    - [fig] 最も相関が高いエッジの経過日数とECvの分布 `resultB_CorrCoef/range_ECv/scat_maxCorrCoef_ECv/`

```
script/B_Visualization_CorrCoef_GeneExp.py
```
- 遺伝子発現値の分布を確認
    - [input] BN result (=ECvデータ) `BayesianNetworkEstimation/CRC/ECv_network_result_0.15.txt`
    - [output] GeneExpノードと腫瘍増殖率の相関係数のリスト `data/data_CorrCoef/GeneExp/`
        - `Parent	Child	CorrelationCoefficients	range_ecv`
    - [fig] 遺伝子発現値のrangeの分布 `resultB_CorrCoef/range_GeneExp/dist_CorrCoef_GeneExp/`



## C. 二値分類予測

1. データ前処理
```
script/C_make_BinaryClassificationMatrix_rangeECv.py
```
- 二値分類予測input matrixを作成
    - [input] BN result (ECv) `BayesianNetworkEstimation/CRC/ECv_network_result_0.15.txt`
    - [input] RawData (GeneExp) `BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt`
    - [input] SelectedEdgesData `data/SelectedEdges_CorrCoef_rangeECv/`
    - [output] 二値分類予測のためのmatrix `data/data_BinaryClassificationDataSet/`

```
script/C_make_BinaryClassificationMatrix_rangeGeneExp.py
```
- 二値分類予測input matrixを作成
    - [input] RawData (GeneExp) `BayesianNetworkEstimation/input_dataset/CRC/CRC_dataset.txt`
    - [input] SelectedEdgesData `data/SelectedEdges_CorrCoef_rangeGeneExpression/`
    - [output] 二値分類予測のためのmatrix `data/data_BinaryClassificationDataSet/`


2. 可視化
```
script/C_BinaryClassification_AllComponent.py
```
- CorrCoef閾値で抽出した全エッジ(orノード)で二値分類予測
    - [fig] ROC curves -> (未完成)
    - [fig] 二値分類予測のためのmatrix `resultC_BinaryClassification/result_AllComponent`

```
script/C_BinaryClassification_AllComponent.py
```
- CorrCoef閾値で抽出した全エッジ(orノード)で二値分類予測
    - [fig] 二値分類予測のためのmatrix `resultC_BinaryClassification/result_SelectedComponent`


## D. Bayesian Network推定 RioEJCA2017 

1. データ前処理
```
script/D_Preprocessing_BNinput_RioEJCA2017.py
```
- BN input matrix作成
    - [input] RawData (RioEJCA2017) `RawData/RioEJCA2017/GSE72970_family.soft`
    - [input] RawData (RioEJCA2017) `RawData/RioEJCA2017/GSE72970_series_matrix.txt`
    - [output] BN input matrix (RioEJCA2017) `BayesianNetworkEstimation/input_dataset/CRC/RioEJCA2017_dataset.txt`
    - [output] 生存時間データ (RioEJCA2017) `data/data_RioEJCA2017/RioEJCA2017_pfsos.txt`


## E. 生存時間解析 RioEJCA2017 

1. データ前処理
```
script/E_make_SurvivalAnalysisMatrix_rangeECv_RioEJCA2017.py
```
- 生存時間解析のためのinput matrixを作成
    - [input] BN result (ECv) `BayesianNetworkEstimation/CRC/ECv_Extrapolation_RioEJCA2019/ECv_extrapolation_RioEJCA2017_dataset_v2.txt`
        - v1 は　前処理済みの生データに過剰にlogをとったもの (ECv range小さい)
        - v2 は　前処理済みデータそのまま
    - [input] 生存時間データ (RioEJCA2017) `data/data_RioEJCA2017/RioEJCA2017_pfsos.txt`
    - [input] SelectedEdgesData `data/SelectedEdges_CorrCoef_rangeECv/`
        - 二値分類の時と同じエッジで生存予測
    - [output] 生存時間解析のためのECvmatrix `data/data_SurvivalAnalysisDataSet_RioEJCA2017/SurvivalAnalysisDataSet_ECv_th06/`
        - multidose : 多剤併用の患者も含む
        - monodose : 単剤患者のみ

2. 可視化
```
script/E_Visualization_Clustering_RioEJCA2017.py
```
- クラスタリング Heatmap
    - [input] 生存時間解析のためのECvmatrix `data/data_SurvivalAnalysisDataSet_RioEJCA2017/SurvivalAnalysisDataSet_ECv_th06/`
    - [output] サンプルごとのクラスターIDとECv `data/data_RioEJCA2017/ClusterIDs/`
    - [fig] クラスタリング結果 `resultE_RioEJCA2017/Dendrogram/`
    - [fig] クラスター数とエッジ数の分布 `resultE_RioEJCA2017/Threshold/`
    - [fig] 階層クラスタリング&Heatmap結果 `resultE_RioEJCA2017/Heatmap/`

```
script/E_Visualization_SurvivalAnalysis_RioEJCA2017.py
```
- クラスタリング Heatmap
    - [input] サンプルごとのクラスターIDとECv `data/data_RioEJCA2017/ClusterIDs/`
    - [input] 生存時間データ `data/data_RioEJCA2017/RioEJCA2017_pfsos.txt`
    - [fig] クラスターごとの生存時間の分布 `resultE_RioEJCA2017/OS/`
    - [fig] 生存時間 `resultE_RioEJCA2017/SurvivalAnalysis/Heatmap/`





## Directory structure
```
.
├── README.md                                              : this file
├── BayesianNetworkEstimation/
│   ├── CRC/
│   ├── SCNEC/
│   └── input_dataset/
│
├── data/                                                  : dataset
│   ├── ArrayDataset/                                      : 
│   ├── data_TumorGrowthRate/
│   ├── data_CorrCoef/
│   ├── data_GeneList/
│   ├── SelectedEdges_CorrCoef_rangeECv/ 
│   ├── SelectedGeneSet_CorrCoef_rangeGeneExpression/
│   ├── data_BinaryClassificationDataSet/
│   ├── data_RioEJCA2017/
│   └── data_SurvivalAnalysisDataSet_RioEJCA2017/
│
├── RawData/                                               : dataset
│   ├── InVivoData/                                        : 
│   └── RioEJCA2017/                                       : 
│
├── resultA_TumorGrowthRate/
│   ├── Clustering/
│   ├── plot_mean/
│   └── plot_median/
│
├── resultB_CorrCoef/
│   ├── range_ECv/
│   └── range_GeneExp/
│
├── resultC_BinaryClassification/
│   ├── result_AllComponent/
│   └── result_SelectedComponent/
│
├── resultE_RioEJCA2017/
│   ├── Dendrogram/
│   ├── Heatmap/
│   ├── Threshold/
│   ├── OS/
│   └── SurvivalAnalysis/
│
└── script/                                                    : script files
    ├── A_data_processing_TumorGrowthRate.py
    ├── B_data_processing.py
    ├── B_deltaECv_cerv54_SCNEC.py
    ├── B_subnetwork_extraction_cerv54_SCNEC.py
    ├── B_Preprocessing_BNinput_CTOS.py
    ├── B_Visualization_CorrCoef_ECv.py
    ├── B_Visualization_CorrCoef_GeneExp.py
    ├── C_make_BinaryClassification_genelist.py
    ├── C_make_BinaryClassificationMatrix_rangeECv.py
    ├── C_make_BinaryClassificationMatrix_rangeGeneExp.py
    ├── C_BinaryClassification_AllComponent.py
    ├── C_BinaryClassification_SelectedComponent.py
    ├── D_Preprocessing_BNinput_RioEJCA2017.py
    ├── E_make_SurvivalAnalysisMatrix_rangeECv_RioEJCA2017.py
    ├── E_Visualization_Clustering_RioEJCA2017.py
    └── E_Visualization_SurvivalAnalysis_RioEJCA2017.py
```



