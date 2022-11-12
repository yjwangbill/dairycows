% Predictions for Kernel Ridge Regression
function Ypred = predictGSAM(Xte, Xtr, kernelFunc, Alp)
  Ktetr = kernelFunc(Xte, Xtr);
  preds = Ktetr * Alp;
  [~, Ypred] = max(preds');
  Ypred = Ypred';