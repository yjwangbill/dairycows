function [alpha, obj] = calcuA_pin(options,K,Y)


%% Initialization
%
if isfield(options,'d')
    d = options.d;
else
    error('Dimension number d not provided.');
end

if isfield(options,'mu')
    mu = options.mu;
else
    mu = 1;
end

if isfield(options,'iter')
    iter = options.iter;
else
    iter = 500;
end

if isfield(options,'lambda')
    lambda = options.lambda;
else
    lambda = 1;
end

if isfield(options,'tau')
    tau = options.tau;
else
    tau = 0;
end
%
n_te = length(Y);
n_tr = size(K,2)/d;

%
%alpha = zeros(n_tr*d, iter);
alpha_cur = zeros(n_tr*d,1)+eps;  %eps  默认值2.2204e-16
gamma = alpha_cur;
obj = [];

%
cL = (n_te*((1+tau)^2/mu)) * max(sum(K.^2,2));
cR = d/mu;
cU = cL+ lambda*cR;  %corresponding Lipschitz constant
%% Calculation
%
for k = 0 : iter
    
    % update delR
    tmpu = 1-Y.*(K*alpha_cur); %1-yf(x)
    u = median([(1+tau)*tmpu/mu, zeros(n_te,1), ones(n_te,1)], 2);  %median（A,2）按行取中位数
    delL = -(1+tau)*K'*(Y.*u)+tau*K'*Y;  %对alpha_cur求导
    
    % update delU
    tmpa = reshape(alpha_cur,n_tr,[]);
    if options.q == 2
        tmpv0 = sqrt(sum(tmpa.^2));
    else
         tmpv0 = max(abs(tmpa));
    end
    tmpv = repmat(mu*max(tmpv0,1),n_tr,1);
    v = alpha_cur ./ reshape(tmpv,[],1);
    delR = v;
    
    % update delF
    delF = delL + lambda*delR;
    
    % update beta： the standard gradient descent solution
    beta = alpha_cur - delF/cU;
    
    if k > 0
        gamma = gamma - 0.5*(k+1)/cU*delF;
    end
    
    % update alpha
    alpha_cur = (2*gamma + (k+1)*beta)/(k+3);
    
    % calculate obj
    if mod(k,50) == 1
%         lossL = sum(tmpu.*u) - 0.5*mu*u'*u;
        lossL = sum(max(tmpu,0));
%         lossR = delR'*v - 0.5*mu*v'*v;
        lossR = sum(tmpv0);
        loss = lossL + lambda*lossR;
        obj = [obj; loss];
    end
    
end


%% Output
alpha = gather(alpha_cur);

end
