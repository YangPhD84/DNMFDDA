# DNMFDDA
To effectively extract low-dimensional key features from high-dimensional, complex data spaces and leverage prior similarity information for predicting potential drug-disease associations, we developed DNMF-DDA, a novel drug repurposing model based on deep non-negative matrix factorization. This model integrates similarity metrics and known associations to extract low-rank representations, enhancing the accuracy of drug-disease association prediction. To improve performance for novel drugs, we employ k-nearest neighbor (KNN) preprocessing to increase the density of prior information within the matrix. Subsequently, two integrated matrices are constructed based on drug and disease similarities, along with optimized association data. During deep matrix factorization, we incorporate graph Laplacian regularization and relaxed constraints to refine local structural features, thereby enhancing the model's capacity to capture complex drug-disease relationships while mitigating the negative impact of insufficient prior knowledge in cold-start scenarios. Additionally, we impose non-negativity constraints to ensure biologically meaningful predictions.

# Requirements
* Matlab >= 2014

# Installation
DNMF-DDA can be downloaded by
```
git clone https://github.com/YangPhD84/DNMFDDA
```
Installation has been tested in a Windows platform.

# Datasets Description
* Gdataset.mat: this file contains information about the gold standard dataset;
* Cdataset.mat: this file contains information about the Cdataset;
* CTDdataset2023.mat: this file contains information about the CTDdataset2023;
* Wdname: disease name;
* Wrname: drug name;
* drug_ChemS: chemical structure similarity matrix;
* drug_AtcS: drug's ATC code similarity matrix;
* drug_SideS: side-effect similarity matrix;
* drug_DDIS: drug-drug interaction similarity matrix;
* drug_TargetS: drug's target profile similarity matrix;
* disease_PhS: disease phenotype similarity matrix;
* disease_DoS: disease ontology similarity matrix;
* didr: disease-drug association matrix.
