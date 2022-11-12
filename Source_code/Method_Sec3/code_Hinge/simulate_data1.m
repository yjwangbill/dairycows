
function [f,X,Y] = simulate_data1(n, p, r, a, b, t,rate) 
    
    U       = a + (b - a).*repmat(rand(n, 1), 1, p);
    W       = a + (b - a).*rand(n, p);
    X       = (W + t*U)/(1 + t);
    f       = @(x_1,x_2) (x_1).^2+(x_2).^2-0.4;
    val     = f(X(:,1),X(:,2));
    s       = sort(abs(val),'ascend');
    label   = find(abs(f(X(:,1),X(:,2)))<=s(0.3*n));
    randm   = randperm(size(label,1));  no_label = ones(n,1);   
    lab_no  = label(randm(1:floor(rate*size(label,1))));
    no_label(lab_no)  =   -1;
    Y       = sign(f(X(:,1),X(:,2))).*(no_label.^r);
end
