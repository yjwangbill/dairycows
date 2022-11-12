function [bestlambda,bestAcc] = GSAM_crossval(X,Y,params)

seqs=2;%交叉验证次数
lambda = params.lambda;
n = length(lambda);
Acc = zeros(n);
bestlambda = 1000000;
bestAcc = 0;

    for j = 1 : n
        params.lambda = lambda(j);
        Acc(j) = crossval(X,Y,seqs,params);              
        if Acc(j) > bestAcc
            bestAcc = Acc(j);
            bestlambda  =  lambda(j);
        end 
    end


