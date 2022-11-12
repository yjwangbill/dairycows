clc,clear
addpath  code_Hinge
n     =   400;               
p     =    40;             
a     =   -1 ;             
b     =    1 ;              
t     =    0.5;         
group = cell(1,p);  
for i = 1:p
    group(1,i)={i};         
end
params.group   = group;
params.order   = 3;         
params.q       = 2;         
params.lambda  = 1000;      
params.tau     = 0.2;       
params.mu      =   1;        
rate           = 0.3;       
num_iter       = 20;
for i = 1:num_iter               
    [~,traindata{i},trainlabel{i}] =  simulate_data(n, p, 1, a, b, t,rate);
    [~,testdata{i},testlabel{i}] = simulate_data(2000, p, 0, a, b, t,rate);
end

%% Traing with Pin loss
for i = 1:num_iter                
         disp(['iteration:  ',num2str(i)])
         [predFunc_Pin, alpha_Pin, fea_Pin,~] = GSAM_pred_Pin(traindata{i}, trainlabel{i}, params);     
         [~, ~,Testing_Pin_Acc(i),~] = GSAM(traindata{i}, trainlabel{i}, testdata{i}, testlabel{i},predFunc_Pin); 
end                   
disp(['Accuracy:  ',num2str(mean(Testing_Pin_Acc))])


