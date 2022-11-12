%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison between GSMAM and others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc
addpath Functions;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n       = 200;                      % sample size of X 
p       = 100;                      % dimension of X                
lambda1 = 1  ;                      % lambda of regularization for threshold
lambda2 =0.01 ;                     % lambda of regularization for mse
r       =   1;                      % sd of noise, s/n ratio = 3
t       =   0;                      % correlation parameter of X
a       = -0.5;                     % lower bound of X
b       = 0.5;                      % upper bound of X
dG      = 1  ;                      % size of group 
nG      = p/dG;                     % number of groups
True_feature=[9,10,99,100];
False_feature=[1:8,11:98,101:p];    
partition=dG*ones(nG,1);partition2=n*ones(nG,1);
cum_part = cumsum(partition);
NIter=25;
tau=ones(nG,1);

Size_fea=zeros(NIter,1);MSE_GSAM=Size_fea;RSSE_GSAM=Size_fea;
TP=zeros(NIter,1);FP=zeros(NIter,1);
CF=zeros(NIter,1);UF=zeros(NIter,1);OF=zeros(NIter,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data generate
for ii=1:NIter
[Xtrain, Ytrain]    = simulate_data(n, p, r, a, b, t);%train data
[Xtest, Ytest]      = simulate_data(n, p, 0, a, b, t);%true/test data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GSMAM method
options.Kernel =  'rbf' ; options.KernelParam=0.5;
start_ind=1;
for i=1:nG
    sel = start_ind:cum_part(i);
    K_train(:,n*(i-1)+1:n*i)=calckernel(options,Xtrain(:,sel),Xtrain(:,sel));
    K_test(:,n*(i-1)+1:n*i)=calckernel(options,Xtrain(:,sel),Xtest(:,sel));
    start_ind = cum_part(i) + 1;
end
para.r=lambda1;
para.kerOpt='Gauss';
para.regOpt='G21';

para.partition=partition2;
[alpha,~,~]=solveMR_linear(K_train,Ytrain,para);
beta=zeros(1,p);
for i=1:p
    beta(i)=norm(alpha(n*(i-1)+1:n*i),2);
end
thre        = sort(beta);  vn=thre(p-4);   %vn : the number of remaining variables
feature     = find(beta > vn);             %Threshold the group sparse vector
K_train2    =[];
K_test2     =[];
%% MSE
for j=1:size(feature,2)
    K_test2(:,n*(j-1)+1:n*j)=calckernel(options,Xtrain(:,feature(j)),Xtest(:,feature(j)));
    K_train2(:,n*(j-1)+1:n*j)=calckernel(options,Xtrain(:,feature(j)),Xtrain(:,feature(j)));
end
partition3=ones(size(feature,2),1);
partition4=size(Xtrain,1)*ones(size(feature,2),1);
para.partition=partition4;
para.r=lambda2;
[alpha2,D,obj]=solveMR_linear(K_train2,Ytrain,para);
ftest=K_test2*alpha2;
MSE(ii)=mean((Ytest-ftest).^2);

%% Show Results
disp(['Niter=',num2str(ii)]);
fprintf('True Feature: %s\n', int2str(True_feature));
fprintf('Selected Feature of GSAM: %s\n', int2str(feature));
Size_fea(ii)=length(feature);
[TP(ii),FP(ii),CF(ii),UF(ii),OF(ii)]=...
    Evalu_Vari_selection(feature,True_feature,False_feature);
    Evalu_Vari_selection(feature,True_feature,False_feature);
end
size_fea=mean(Size_fea);tp=mean(TP);fp=mean(FP);
cf=sum(CF);uf=sum(UF);of=sum(OF);
mse=mean(MSE);
disp(['Size=',num2str(size_fea),' TP=',num2str(tp), ' FP=',num2str(fp)]);  
disp(['C=',num2str(cf), ' U=',num2str(uf), ' O=',num2str(of), ' mse=',num2str(mse)]); 