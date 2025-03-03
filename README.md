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

# Functions Description
* approx_seminmf: Matrix initialization code based on BasicNMF(*) function;
* code: The main function code of the DNMF-DDA model;
* Result: Disease scores and OMIM codes for the top 10 predictive scores for all drug associations were recorded; 
* DNMFDDA_Demo.m: DNMF-DDA model demonstration code;

# Functions Description
* approx_seminmf: Matrix initialization code based on BasicNMF(*) function;
* code: The main function code of the DNMF-DDA model;
* Result: Disease scores and OMIM codes for the top 10 predictive scores for all drug associations were recorded; 
* DNMFDDA_Demo.m: DNMF-DDA model demonstration code;
* KNN_diseaseS.m: This function fills the association matrix with all zero rows using disease similarity;
* KNN_drugS.m: This function fills the association matrix with all zero columns using drug similarity;
* fDRNMF.m: This function is the main function of DNMF-DDA;
* Wupdate.m: This function is used to update W in DNMF-DDA;
* Zupdate.m: This function is used to update Z in DNMF-DDA;
* fdrinit.m: This function is used to initialize the base matrix and coefficient matrix in deep matrix decomposition;

# Instructions
We provide detailed step-by-step instructions for running DNMF-DDA model.

**Step 1**: add datasets\functions paths
```
addpath('Datasets');
addpath('approx_seminmf');
addpath('code');
load Gdataset              
% load Cdataset             
% load CTDdataset2023

%% Selecting dataset
vars = whos;
T = 0;
for i = 1:length(vars)
    if isequal(vars(i).size(1), 313)
        T = 1;
        break;
    elseif isequal(vars(i).size(1), 409)
        T = 2;
        break;
    elseif isequal(vars(i).size(1), 278)
        T = 3;
        break;
    end
end
```
**Step 2**: load association matirx and similarity matrices

```
drug_ChemS=(drug_ChemS+drug_ChemS')/2;
drug_AtcS=(drug_AtcS+drug_AtcS')/2;
drug_SideS=(drug_SideS+drug_SideS')/2;
drug_DDIS=(drug_DDIS+drug_DDIS')/2;
drug_TargetS=(drug_TargetS+drug_TargetS')/2;
disease_PhS=(disease_PhS+disease_PhS')/2;
disease_DoS=(disease_DoS+disease_DoS')/2;

Wrr = (drug_ChemS + drug_AtcS + drug_SideS + drug_DDIS + drug_TargetS)/5;
Wdd = (disease_PhS + disease_DoS)/2;

Wdr=didr;
Wrd = Wdr';
[dn,dr] = size(Wdr);
```
**Step 3**: parameter Settings

The hyper-parameters are fixed.
```
MaxIter = 400;      
tol1 = 5*1e-4;      
tol2 = 5*1e-7;      
alpha = 0.01;
beta = 1;
lammad = 1;
gamma = 1;
```
**Step 4**: KNN preprocessing of all-zero rows or all-zero columns of the association matrix
```
K = 10;
row_no=find(sum(Wdr,2)==0);
if isempty(row_no)==0
     P_TMat_new0= KNN_diseaseS(Wdr,Wdd,K);
     P_TMat_dd= P_TMat_new0+Wdr;
else
    P_TMat_dd=Wdr;
end
Drug_total_Association = [Wrr,P_TMat_dd'];
field_Drug = {Drug_total_Association};
Drug_similarity = {Wrr};

col_no=find(sum(Wdr,1)==0);
if isempty(col_no)==0
    P_TMat_new00= KNN_drugS(Wdr,Wrr,K);
    P_TMat_rr= P_TMat_new00+Wdr;
else
    P_TMat_rr=Wdr;
end

Disease_total_Association = [Wdd,P_TMat_rr];
field_Disease = {Disease_total_Association};
Disease_similarity = {Wdd};
```
**Step 4**: Set the number and size of decomposition layers and execute the fDRDMF function.
```
[R_m,R_n] = size(field_Drug{1,1});
R_min_mn = min(R_m,R_n);
R_layer_2 = floor(R_min_mn*0.8);
R_layer_4 = floor(R_min_mn*0.6);
R_layers = [R_layer_2 R_layer_4];

[D_m,D_n] = size(field_Disease{1,1});
D_min_mn = min(D_m,D_n);
D_layer_2 = floor(D_min_mn*0.8);
D_layer_4 = floor(D_min_mn*0.6);

D_layers = [D_layer_2 D_layer_4];

[R_M, R_Z, R_W] = fDRDMF(Wdr,field_Drug,Drug_similarity,R_layers,alpha,beta,lammad,gamma,MaxIter,tol1,tol2,1,T);

[D_M, D_Z, D_W] = fDRDMF(Wdr,field_Disease,Disease_similarity,D_layers,alpha,beta,lammad,gamma,MaxIter,tol1,tol2,0,T);
```
# A Quickstart Guide
Users can immediately start playing with DNMF-DDA running ``` DNMFDDA_Demo.m ``` in matlab.
* ```DNMFDDA_Demo.m```: it demonstrates a process of predicting drug-disease associations on the gold standard dataset (Gdataset) by DNMF-DDA algorithm.

# Run DNMF-DDA on User's Own Data
We provided instructions on implementing DNMF-DDA model with user's own data. One could directly run DNMF-DDA model in ``` DNMFDDA_Demo.m```  with custom data by the following instructions.

**Step 1**: Prepare your own data and add the corresponding dataset files.

The required data includes drug-disease association matirx and similarity matrices, which are all saved by ```mat``` files.

**Step 2**: Modify the drug and disease information in ```DNMFDDA_Demo.m```

# Contact
If you have any questions or suggestions with the code, please let us know. Contact Mengyun Yang at mengyun_yang@126.com
