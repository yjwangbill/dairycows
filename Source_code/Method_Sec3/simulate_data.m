
function [f,X,Y] = simulate_data(n, p, r, a, b, t,rate) 
    
    U       = a + (b - a).*repmat(rand(n, 1), 1, p);
    W       = a + (b - a).*rand(n, p);
    X       = (W + t*U)/(1 + t);
    f       = @(x_1,x_2,x_3,x_4,x_5,x_6) x_1+x_2.^2+x_3.^3+sin(pi*x_4)-log(x_5+3)-abs(x_6)+1;
    val     = f(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6));
    s       = sort(abs(val),'ascend');
    label   = find(abs(f(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6)))<=s(0.3*n));
    randm   = randperm(size(label,1));  no_label = ones(n,1);   
    lab_no  = label(randm(1:floor(rate*size(label,1))));
    no_label(lab_no)  =   -1;
    Y       = sign(f(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6))).*(no_label.^r);
end
