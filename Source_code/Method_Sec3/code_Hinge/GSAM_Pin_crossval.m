function [best_tau,bestlambda,bestAcc] = GSAM_Pin_crossval(X,Y,params)

seqs   = params.numPartsKFoldCV;%交叉验证次数
lambda = params.lambda;
tau    = params.tau;
n = length(lambda);
m = length(tau);
Acc = zeros(n);
bestlambda = 1000000;
bestAcc = 0;
   for i = 1 : m
       params.tau   =   tau(i);
    for j = 1 : n
        params.lambda = lambda(j);
        Acc(j,i) = crossval_pin(X,Y,seqs,params);              
    end
   end
[bestAcc,II] = max(max(Acc));
[~,I]        = max(Acc(:,II));
bestlambda   = lambda(I);
best_tau     = tau(II);

