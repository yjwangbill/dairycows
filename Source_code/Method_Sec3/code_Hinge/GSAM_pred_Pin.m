function [predFunc,kernelFunc, alpha, fea] = GSAM_pred_Pin(X, Y, params)

% this function impletes GSAM_pred : Group Sparse Shrunk Additive Model
% input: X trainData n*D
% 	     Y trainLabel n*1: +1 for positive class, -1 for negative class
%
% params: a structure which optionally contain the following parameters
% 		- numPartsKFoldCV
% 		- numTrailKFoldCV
% 		- numLambdaCands
% 		- lambdaRange
%       - lambda : use certain lambda or use direct to tune hyperpara
% Alp the coefficient vector (Alp)
%
% output: predFunc: A function handle which can be used to estimate the function for new points
% 		  

KernelOpts.setting = 'espKernel';
if isfield(params,'bws')
    KernelOpts.bws = params.bws;
end
if isfield(params,'group')
    KernelOpts.group = params.group;
end
if isfield(params,'order')
    KernelOpts.order = params.order;
else
    KernelOpts.order = 2;
end

kernelFunc = kernelSetup(X, Y, KernelOpts);
K = kernelFunc(X, X);

unic = unique(Y); %returns the same values as in Y but with no repetitions. unic will be sorted. 
c = length(unic); % 二分类c=2;
n = length(Y);  % n应该是样本个数
d = size(K, 2)/n;  %核矩阵大小的一半，样本个数的一半，变量选择后的维数
alpha = zeros(d*n, c);

options.d = d;
options.lambda = params.lambda;
options.tau    =  params.tau;
options.mu     = params.mu;
options.q      = params.q;
for t = 1 : c
    Yt = -ones(n, 1);
    Yt(Y == unic(t)) = 1;  %分别找到标签为-1/1的数据
    [alpha(:,t)] = calcuA_pin(options, K, Yt);  %calcuA ????
    tmp = reshape(alpha(:,t),n,[]);
    tmpfea(t,:) = sqrt(sum(tmp.^2));
end

fea = tmpfea(1,:);

predFunc = @(Xte, Yte) predictGSAM(Xte, X, kernelFunc, alpha);

end

% Predictions for Kernel Ridge Regression
function Ypred = predictGSAM(Xte, Xtr, kernelFunc, Alp)
  Ktetr = kernelFunc(Xte, Xtr);
  preds = Ktetr * Alp;
  [~, Ypred] = max(preds');
  Ypred = Ypred';
end