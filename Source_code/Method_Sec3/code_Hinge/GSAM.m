function [F,pre, Traning_Acc, Testing_Acc,recall,TrYpred, TeYpred] = GSAM(traindata, trainlabel, testdata, testlabel,predFunc)
TrYpred = predFunc(traindata, trainlabel);
TrYpred(TrYpred==1) = -1;
TrYpred(TrYpred==2) = 1;
Traning_Acc = sum(TrYpred==trainlabel)/size(trainlabel,1);

TeYpred = predFunc(testdata, testlabel);
TeYpred(TeYpred==1) =  -1;
TeYpred(TeYpred==2) =   1;
Testing_Acc = sum(TeYpred==testlabel)/size(testlabel,1);

index_tp = find(TeYpred==1);
index_fp = find(TeYpred == -1);
TP = sum(TeYpred(index_tp)==testlabel(index_tp));
FP = sum(TeYpred(index_tp)~=testlabel(index_tp));
TN = sum(TeYpred(index_fp)==testlabel(index_fp));
FN = sum(TeYpred(index_fp)~=testlabel(index_fp));
recall      = TP/(TP+FN);
pre        = TP/(TP+FP);
F = 2*pre*recall/(pre+recall);


