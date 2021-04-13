# CTOS

## Directory structure
```
.
├── README.md                                              : this file
├── BayesianNetworkEstimation/
│   ├── CRC/
│   ├── SCNEC/
│   └── input_dataset/
├── data/                                                  : dataset
│   ├── ArrayDataset/                                      : 
│   ├── InVivoData/                                        : 
│   ├── SelectedEdges_CorrCoef_rangeECv/ 
│   └── SelectedGeneSet_CorrCoef_rangeGeneExpression/
├── data_BinaryClassification/
├── data_CRC/
├── data_TGR/
│   └── TGR.csv
├── resultA_CorrCoef/
│   ├── range_ECv/
│   └── range_GeneExp/
├── resultA_TumorGrowthRate/
│   ├── Clustering_TGR/
│   ├── plot_mean/
│   └── plot_median/
├── resultB_BinaryClassification/
│   ├── Classification_result/
│   └── ROC_curves/
├── resultC_BinaryClassification_SelectedComponent/
│   └── Classification_result/
└── script/                                                    : script files
    ├── Preprocessing_FeatureExtraction_CRC.py
    ├── Preprocessing_Labeling_TGR_CRC.py
    ├── Visualization_Clustering_TGR_CRC.py
    ├── Visualization_ECv.py
    ├── Visualization_GeneExp.py
    ├── Preprocessing_binaryclassification_all.py
    ├── BinaryClassification_all.py
    ├── Preprocessing_binaryclassification_SelectedComponent.py
    ├── BinaryClassification_SelectedComponent.py
    ├── data_processing.py
    ├── deltaECv_cerv54_SCNEC.py
    └── subnetwork_extraction_cerv54_SCNEC.py
```
