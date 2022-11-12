
function Acc0 = crossval(X,Y,seqs,params)
%seqs:������֤����
[kx,n] = size(X);
[ky,m] = size(Y);
if kx ~= ky, disp('Incompatible X and Y'); return; 
else k = kx; end

Xorig = X;
Yorig = Y;
seqlen = ceil(k/seqs);%ceil ��������ȡ��
for i = 1:seqs
   begin = (i-1)*seqlen+1;
   endin = min(i*seqlen,k);
   Xtrain{i} = Xorig([1:begin-1,endin+1:k],:);%ȡ������֤��ѵ������
   Ytrain{i} = Yorig([1:begin-1,endin+1:k],:);
   Xtest{i} = Xorig(begin:endin,:);%ȡ������֤�Ĳ�������
   Ytest{i} = Yorig(begin:endin,:);
   
   [predFunc_ATk, kernelFunc_ATk,alpha_ATk, fea_ATk] = GSAM_pred(Xtrain{i}, Ytrain{i}, params);
   Ypred_ATk{i} = predFunc_ATk(Xtest{i}, Ytest{i});
   Ypred_ATk{i}(find(Ypred_ATk{i}==1)) = -1;
   Ypred_ATk{i}(find(Ypred_ATk{i}==2)) =1;
   Testing_ATk_Acc(i) = sum(Ypred_ATk{i}==Ytest{i})/size(Ytest{i},1);%������֤�����
   % [trainPredictY,testPredictY,traErr(i),TstErr(i)]=ComputeKernel(Xtrain,Ytrain,Xtest,Ytest,rbf_var,lamda);%������֤�����
end
Acc0=mean(Testing_ATk_Acc);
