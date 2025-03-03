clc,clear all

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

%%
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

MaxIter = 400;      
tol1 = 5*1e-4;     
tol2 = 5*1e-7;      
alpha = 0.01;
beta = 1;
lammad = 1;
gamma = 1;

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

%% DRDMF algorithm
%drug_DRDMF
[R_m,R_n] = size(field_Drug{1,1});
R_min_mn = min(R_m,R_n);
R_layer_2 = floor(R_min_mn*0.8);
R_layer_4 = floor(R_min_mn*0.6);
R_layers = [R_layer_2 R_layer_4];

%disease_DRDMF
[D_m,D_n] = size(field_Disease{1,1});
D_min_mn = min(D_m,D_n);
D_layer_2 = floor(D_min_mn*0.8);
D_layer_4 = floor(D_min_mn*0.6);

D_layers = [D_layer_2 D_layer_4];

[R_M, R_Z, R_W] = fDRDMF(Wdr,field_Drug,Drug_similarity,R_layers,alpha,beta,lammad,gamma,MaxIter,tol1,tol2,1,T);

[D_M, D_Z, D_W] = fDRDMF(Wdr,field_Disease,Disease_similarity,D_layers,alpha,beta,lammad,gamma,MaxIter,tol1,tol2,0,T);

if T == 1
    M_Result = (R_M(:,594:end)' + D_M(:,314:end))/2 ;% Gdataset
elseif T == 2
    M_ResultMat = (R_M(:,664:end)' + D_M(:,410:end))/2 ;% Cdataset
elseif T == 3
    M_ResultMat = (R_M(:,1238:end)' + D_M(:,279:end))/2 ;% CTDdataset
end
