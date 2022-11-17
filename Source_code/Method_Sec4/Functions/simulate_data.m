
function [X,Y] = simulate_data(n, p, r, a, b, t)
    U       = a + (b - a).*repmat(rand(n, 1), 1, p);
    W       = a + (b - a).*rand(n, p);
    X       = (W + t*U)/(1 + t);
    Y       = true_function(X) +r*random('chi',1,n,1);
end

%Example 1--Function 1
function Y = true_function(X)

  Y=-2*sin(2.*X(:,9)) + 3*X(:,10).^2 + 2*sin(X(:,11))./(2-sin(X(:,11)))-exp(-X(:,12))+X(:,97).^3+1.5*(X(:,97)-1).^2+X(:,98)+5.*sin(exp(-0.5.*X(:,99)))-5.*normcdf(X(:,100),0.5,0.8);
end
