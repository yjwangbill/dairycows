function [W,D, obj] = solveMR_linear(X, Y, para)

% min_{w} 1/(n*sigma) * \sum_{i=1}^n (1 - K((y_i - W'*x_i)/sigma))
% solve the objecitive in modal regression via half quadratic optimization
% min_{w, D} ||\sqrt{D}(Y-X'*w)||_2^2 + \sum_{i=1}^n \psi(D(i,i))
%
% choice of kernel in para.kerOpt
%       - Epanech: Epanechnikov kernel
%       - Gauss: Gaussian kernel
%       - Quartic
%       - Triangular
%       - Logistic
%       - Sigmoid
%
% choice of regularization term of W in para.reg
%       - Fro: ||W||_F^2
%       - G21: ||W||_{G2,1}
%       - L1: ||W||_1
%
% Author: Xiaoqian Wang
%

if isfield(para,'maxIter')
    maxIter = para.maxIter;
else
    maxIter = 10;
end

if isfield(para,'r')
    lambda = para.r;
else
    lambda =0.1;
end

if isfield(para,'kerOpt')
    kerOpt = para.kerOpt;
else
    kerOpt = 'Gauss';
end

if isfield(para,'regOpt')
    regOpt = para.regOpt;
else
    regOpt = 'Fro';
end

if isfield(para,'partition')
    partition = para.partition;
    %regOpt = 'G21';
else 
    partition = [];
end

[n, d] = size(X);
c = size(Y,2);
kerSet = {'Gauss','Logistic','Sigmoid'};
sigOpt = sum(strcmp(kerOpt,kerSet));

% initialize W
% W = zeros(d, c);
W = rand(d, c);

% initialize sigma
diff = Y-X*W;
if sigOpt == 1
%sigma =50;
%sigma=2.5;
sigma =  sqrt(sum(diff.*diff)/(0.2*n));
else
    sigma = max(abs(diff));
end
% sigma = set_sigma(diff);

for it = 1 : maxIter
    
    % update D
    cmdD = ['D = comD_', kerOpt, '(sigma,diff);'];
    eval(cmdD);
    D = D+eps;
    
    % update w
    cmdW = ['[W, regV] = comW_', regOpt, '(X,Y,D,lambda,partition);'];
    eval(cmdW)
    
    % update sigma
    diff = Y-X*W;
%     sigma = set_sigma(diff);
%     sigma = 2*(sum(diff.*diff)./(n));

    if sigOpt == 1
 % sigma =2.5;
  % sigma=1.5;
sigma =  sqrt(sum(diff.*diff)/(0.2*n));
    else
        sigma = max(abs(diff));
    end

    % update obj
    tmp = diff./repmat(sigma,n,1);
    cmdObj = ['tmpobj = comObj_', kerOpt, '(tmp);'];
    eval(cmdObj);
   obj(it) = n./sigma-(tmpobj./sigma - lambda*regV);   %%argmin
  % obj(it) = tmpobj./sigma - lambda*regV;               %%argmax
%     fprintf('%d\n', obj(it));
    
    if it>1
        if abs(obj(it-1) - obj(it))<1e-5*obj(it-1)
            break;
        end
    end
end

end

%%% functions for computing D with different choice of kernel %%%
function D = comD_Epanech(tot_sigma,tot_diff)
    for ii = 1 : size(tot_diff,2)
        diff = tot_diff(:,ii);
        sigma = tot_sigma(ii);
        idx = (abs(diff/sigma) <= 1);
        D(:,ii) = 0.75/(sigma^1)*idx;
    end
end

function D = comD_Quartic(tot_sigma,tot_diff)
    for ii = 1 : size(tot_diff,2)
        diff = tot_diff(:,ii);
        sigma = tot_sigma(ii);
        tmp = diff/sigma;
        idx = (abs(tmp) <= 1);
        D(:,ii) = 15/(8*sigma^1)*(1-tmp.^2).*idx;
    end
end

function D = comD_Triangular(tot_sigma,tot_diff)
    for ii = 1 : size(tot_diff,2)
        diff = tot_diff(:,ii);
        sigma = tot_sigma(ii);
        tmp = diff/sigma;
        idx = (abs(tmp) <= 1);
        D(:,ii) = 0.5/(sigma^1).*idx./(eps+abs(tmp));
    end
end

function D = comD_Gauss(tot_sigma,tot_diff)
    for ii = 1 : size(tot_diff,2)
        diff = tot_diff(:,ii);
        sigma = tot_sigma(ii);
        tmp = diff/sigma;
        D(:,ii) = 0.5/(sigma^1)*exp(-0.5*(diff.*diff)/(sigma*sigma));
    end
end

function D = comD_Logistic(tot_sigma,tot_diff)
    for ii = 1 : size(tot_diff,2)
        diff = tot_diff(:,ii);
        sigma = tot_sigma(ii);
        tmp = diff/sigma;
        D(:,ii) = 0.5./(sigma^1*tmp+eps).*(-exp(-tmp)+exp(tmp))./...
        ((2+exp(-tmp)+exp(tmp)).^2);
    end
end

function D = comD_Sigmoid(tot_sigma,tot_diff)
    for ii = 1 : size(tot_diff,2)
        diff = tot_diff(:,ii);
        sigma = tot_sigma(ii);
        tmp = diff/sigma;
        D(:,ii) = 1./(pi*sigma^1*tmp+eps).*(-exp(-tmp)+exp(tmp))./...
        ((exp(-tmp)+exp(tmp)).^2);
    end
end


%%% functions for computing W with different choice of regularization %%%
function [W, regValue] = comW_Fro(X,Y,totD,lambda,partition)
    d = size(X,2);
    for ii = 1 : size(totD,2)
        D = diag(totD(:,ii));
        W(:,ii) = pinv(X'*D*X+lambda*eye(d))*X'*D*Y(:,ii);
    end
    regValue = norm(W,'fro')^2;
end

function [W, regValue] = comW_L1(X,Y,D,lambda,partition)
    d = size(X,2);
    b = D.*Y;
    for ii = 1 : size(D,2)
        A = repmat(D(:,ii),1,d).*X;
        W(:,ii) = feature_sign(A, b(:,ii), lambda);
    end
    regValue = sum(sum(abs(W)));
end

function [W, regValue] = comW_G21(X,Y,D,lambda,partition)
    d = size(X,2);
    A = [];
    for ii = 1 : size(D,2)
        A = blkdiag(A,repmat(D(:,ii),1,d).*X);
    end
    b = reshape(D.*Y,[],1);
    tau = ones(size(partition,1),1);
    [z, regValue] = Gene_BSR(A, b, lambda, partition, tau);
    W = reshape(z,d,[]);
end

function [W, regValue] = comW2_G21(X,Y,D,lambda,partition)
    d = size(X,2);
    tau = ones(size(partition,1),1);
    b = D.*Y;
    for ii = 1 : size(D,2)
        A = repmat(D(:,ii),1,d).*X;
        [W(:,ii), reg(ii)] = Gene_BSR(A, b(:,ii), lambda, partition, tau);
    end
    regValue = sum(reg);
end


%%% functions for computing objective with different choice of kernel %%%
function obj = comObj_Epanech(u)
    obj = sum(max(0,0.75*(1-u.^2)));
end

function obj = comObj_Quartic(u)
    obj = 15/16*sum(max(0,1-u.^2).^2);
end

function obj = comObj_Triangular(u)
    obj = sum(max(0,1-abs(u)));
end

function obj = comObj_Gauss(u)
    obj = sum(exp(-0.5*u.^2));
end

function obj = comObj_Logistic(u)
    obj = sum(1./(2+exp(-u)+exp(u)));
end

function obj = comObj_Sigmoid(u)
    obj = sum(2./(pi*(exp(-u)+exp(u))));
end