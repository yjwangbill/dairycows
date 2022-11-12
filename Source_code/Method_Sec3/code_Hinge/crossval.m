
function Acc0 = crossval(X,Y,seqs,params)
%seqs:交叉验证次数
[kx,n] = size(X);
[ky,m] = size(Y);
if kx ~= ky, disp('Incompatible X and Y'); return; 
else k = kx; end

Xorig = X;
Yorig = Y;
seqlen = ceil(k/seqs);%ceil 向正方向取整
for i = 1:seqs
   begin = (i-1)*seqlen+1;
   endin = min(i*seqlen,k);
   Xtrain{i} = Xorig([1:begin-1,endin+1:k],:);%取交叉验证的训练样本
   Ytrain{i} = Yorig([1:begin-1,endin+1:k],:);
   Xtest{i} = Xorig(begin:endin,:);%取交叉验证的测试样本
   Ytest{i} = Yorig(begin:endin,:);
   
   [predFunc_ATk, kernelFunc_ATk,alpha_ATk, fea_ATk] = GSAM_pred(Xtrain{i}, Ytrain{i}, params);
   Ypred_ATk{i} = predFunc_ATk(Xtest{i}, Ytest{i});
   Ypred_ATk{i}(find(Ypred_ATk{i}==1)) = -1;
   Ypred_ATk{i}(find(Ypred_ATk{i}==2)) =1;
   Testing_ATk_Acc(i) = sum(Ypred_ATk{i}==Ytest{i})/size(Ytest{i},1);%交叉验证的误差
   % [trainPredictY,testPredictY,traErr(i),TstErr(i)]=ComputeKernel(Xtrain,Ytrain,Xtest,Ytest,rbf_var,lamda);%交叉验证的误差
end
Acc0=mean(Testing_ATk_Acc);
