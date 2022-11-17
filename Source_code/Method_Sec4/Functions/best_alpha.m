function [feature,mse_test]=best_alpha(Xtrain,Ytrain,Xtest,Ytest,para,n,p,sign,option)
sel=1;
K_train    =[];
K_test     =[];
options = option;
for i=1:p
    K_train(:,n*(i-1)+1:n*i)=calckernel(options,Xtrain(:,sel),Xtrain(:,sel));
    K_test(:,n*(i-1)+1:n*i)=calckernel(options,Xtrain(:,sel),Xtest(:,sel));
    sel =sel + 1;
end

if sign==2
      para.partition=n*ones(p,1);
end
[alpha,~,~]=solveMR_linear(K_train,Ytrain,para);
beta=zeros(1,p);
for i=1:p
    beta(i)=norm(alpha(n*(i-1)+1:n*i),sign);
end
thresh=sort(beta);   vn=thresh(p-8);      %  vn : the number of remaining variables
feature   = find(beta > vn);              %  Threshold the sparse vector by vn
if isempty(feature)
feature=[1];
end
disp(['NIter_lambda=',num2str(para.r)]);

%% testing
if sign == 2
   partition_valid=size(Xtrain,1)*ones(size(feature,2),1);
   para.partition=partition_valid;
end
ftest=K_test*alpha;
mse_test=mean((Ytest-ftest).^2);

