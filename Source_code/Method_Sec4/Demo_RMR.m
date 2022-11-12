%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison between GSMAM and others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc
addpath Functions;
%parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n       = 200;                              % sample size of X 
p       = 100;                              % dimension of X                               
r       = 1  ;                              % sd of noise, s/n ratio = 3
t       =   0;                              % correlation parameter of X
a       = -0.5;                             % lower bound of X
b       = 0.5;                              % upper bound of X

True_feature=[9,10,99,100];
False_feature=[1:8,11:98,101:p];      
NIter=10;

Size_fea=zeros(NIter,1);
TP=zeros(NIter,1);FP=zeros(NIter,1);
CF=zeros(NIter,1);UF=zeros(NIter,1);OF=zeros(NIter,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data generate

for ii=1:NIter
[Xtrain, Ytrain]    = simulate_data(n, p, r, a, b, t);%train data
[Xtest, Ytest]      = simulate_data(n, p, 0, a, b, t);%true/test data
Xtrain = scaleData(Xtrain);
Xtest = scaleData(Xtest);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GSMAM method
para.r   =  0.5;                 % parameter of regularization;
para.kerOpt='Gauss';             % modal kernel
para.regOpt='L1';                
[alpha,~,~]=solveMR_linear(Xtrain,Ytrain,para);
alpha1=alpha.*(abs(alpha)>=0.5);   %Threshold the group sparse vector
feature     = find( alpha1 > eps);
feature     = feature';
fprintf('True Feature: %s\n', int2str(True_feature));
fprintf('Selected Feature of RMR: %s\n', int2str(feature));
%% mse caculation
Xtest2=Xtest(:,feature);
Xtrain2=Xtrain(:,feature);
[alpha2,D,obj]=solveMR_linear(Xtrain2,Ytrain,para);
 ftest=Xtest2*alpha2;  
%% Show Results
Size_fea(ii)=length(feature);
[TP(ii),FP(ii),CF(ii),UF(ii),OF(ii)]=...
Evalu_Vari_selection(feature,True_feature,False_feature);
MSE(ii)=sum((Ytest-ftest).^2)/size(Ytest,1);  

end
size_fea=mean(Size_fea);tp=mean(TP);fp=mean(FP);
cf=sum(CF);uf=sum(UF);of=sum(OF);
mse=mean(MSE);


disp(['MSE=',num2str(mse)]);  
disp(['C=',num2str(cf), ' U=',num2str(uf), ' O=',num2str(of)]); 
disp(['size=',num2str(size_fea), ' TP=',num2str(tp), ' FP=',num2str(fp)]); 
