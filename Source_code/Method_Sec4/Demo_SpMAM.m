%%---------------------------------------------------------- 
% comparison between SpMAM and others
%%---------------------------------------------------------- 
close all;clear;clc
addpath Functions;

%%---------------------------------------------------------- 
%% parameter setting
NIter            =        100    ;           
n                =        200    ;                      % sample size of X 
p                =        100    ;                      % dimension of X            
r                =        0.5    ;                      % sd of noise, s/n ratio = 3
lambda           =        0.2    ;                      % lambda of regularization for threshold
t                =         0     ;                      % correlation parameter of X
a                =        -0.5   ;                      % lower bound of X
b                =         0.5   ;                      % upper bound of X
norm             =         2     ;               
True_feature     =  [9,10,11,12,97,98,99,100];
False_feature    =  [1:8,13:96,101:p];      
Size_fea         =  zeros(NIter,1);

TP      =   zeros(NIter,1); 
FP      =   zeros(NIter,1);
CF      =   zeros(NIter,1);
UF      =   zeros(NIter,1);
OF      =   zeros(NIter,1);

options.Kernel       =  'rbf';  
options.KernelParam  =   0.5 ;                          %bandwidth of RKHS

para.r       =    lambda;
para.kerOpt  =   'Gauss';                               %modal kernel
para.regOpt  =      'L1';                  

%%%%---------------------------------------------------------- 
%data generate
for ii=1:NIter
[Xtrain, Ytrain]    =    simulate_data(n, p, r, a, b, t);    %train data
[Xtest,   Ytest]    =    simulate_data(n, p, 0, a, b, t);    %true/test data
%%---------------------------------------------------------- 
%SpMAM method
%% train and test
[feature,MSE(ii)] = best_alpha(Xtrain,Ytrain,Xtest,Ytest,para,n,p,norm,options);
 %% Show Results
disp(['Niter=',num2str(ii)]);
fprintf('True Feature: %s\n', int2str(True_feature));
fprintf('Selected Feature of GSAM: %s\n', int2str(feature));
Size_fea(ii)=length(feature);
[TP(ii),FP(ii),CF(ii),UF(ii),OF(ii)]=Evalu_Vari_selection(feature,True_feature,False_feature);
end
size_fea=mean(Size_fea);tp=mean(TP);fp=mean(FP);
cf=sum(CF);uf=sum(UF);of=sum(OF);
ASE=mean(MSE);
