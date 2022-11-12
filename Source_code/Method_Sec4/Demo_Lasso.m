%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo_Lasso
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc
addpath Data Functions;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters
n       = 200;                      % sample size of X 
p       = 100;                      % dimension of X
r       =   1;                      % sd of noise, s/n ratio = 3                
t       =   0;                      % correlation parameter of X
a       = -0.5;                     % lower bound of X
b       = 0.5;                      % upper bound of X
True_feature=[9:1:10,99:1:100];
False_feature=[1:8, 11:98,101:p];      
NIter=10;                          %iteration
lambda=1;                          %parameter of regularization  

Size_fea=zeros(NIter,1);MSE_GSAM=Size_fea;RSSE_GSAM=Size_fea;
TP=zeros(NIter,1);FP=zeros(NIter,1);
CF=zeros(NIter,1);UF=zeros(NIter,1);OF=zeros(NIter,1);

for ii=1:NIter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation data
% training data set
[Xtrain, Ytrain]    = simulate_data(n, p, r, a, b, t);
% test data set
[Xtest, Ytest]      = simulate_data(n, p, 0, a, b, t);
Xtrain = scaleData(Xtrain);      %scale
Xtest = scaleData(Xtest);        %scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Method : Lasso
alpha=feature_sign(Xtrain,Ytrain,lambda);
feature = find(abs(alpha)>=1);
feature=feature';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xtest2=Xtest(:,feature);
Xtrain2=Xtrain(:,feature);
alpha2=feature_sign(Xtrain2,Ytrain,lambda);
ftest=Xtest2*alpha2;  
%% Show Results
fprintf('True Feature: %s\n', int2str(True_feature));
fprintf('Selected Feature of Lasso: %s\n', int2str(feature));
Size_fea(ii)=length(feature);
[TP(ii),FP(ii),CF(ii),UF(ii),OF(ii)]=...
    Evalu_Vari_selection(feature,True_feature,False_feature);
  MSE(ii)=sum((Ytest-ftest).^2)/size(Ytest,1);
end
size_fea=mean(Size_fea);tp=mean(TP);fp=mean(FP);
cf=sum(CF);uf=sum(UF);of=sum(OF);
mse=mean(MSE);
disp(['Size=',num2str(size_fea),' TP=',num2str(tp), ' FP=',num2str(fp)]);  
disp(['C=',num2str(cf), ' U=',num2str(uf), ' O=',num2str(of), ' mse=',num2str(mse)]); 
