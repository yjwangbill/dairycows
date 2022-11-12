
function [X,Y] = simulate_data(n, p, r, a, b, t)
    U       = a + (b - a).*repmat(rand(n, 1), 1, p);
    W       = a + (b - a).*rand(n, p);
    X       = (W + t*U)/(1 + t);
    Y       = true_function(X) +r*random('chi',1,n,1);
    %r*normrnd(0,1,n,1);
    %r*random('chi',1,n,1);
    %random('chi',2,n,1);exprnd(n,1)
    %noise=exprnd(1,size(tr_y0));
    %noise=r*nctrnd(3,0,n,1);
end

%Example 1--Function 1
function Y = true_function(X)

  Y=-2*sin(2.*X(:,9)) + 5*X(:,10).^2-1/3 + 2*sin(X(:,11))./(2-sin(X(:,11)))-0.5+exp(-X(:,12))+X(:,97).^3+1.5*(X(:,97)-1).^2+X(:,98)+5.*sin(exp(-0.5.*X(:,99)))-5.*normcdf(X(:,100),0.5,0.8);
end

% %Example 2--Function 2
% function Y = true_function(X)
%    %Y=-2*sin(2*X(:,9))+5*X(:,10).^2-1/3 + sin(X(:,99))-0.5 + exp(-X(:,100))+exp(-1)-1;
% end
