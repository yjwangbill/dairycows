function K = espKernels(X, Y, bws)
% Computes Kernels using elementary symmetric polynomials.
% X :testMatrix, Y : training matrix
% return K num_te*(num_tr * (dimension-1)dimension/2 )
% Author: Xiaoqian Wang

  % prelims
  [num_tr, dimension] = size(Y);
  num_te = size(X, 1);
  order = 2;
  
  % Now construct the ESP kernels
  com = nchoosek(1:dimension,order);
  Y1 = reshape(Y(:,com(:,1)),1,[]);
  Y2 = reshape(Y(:,com(:,2)),1,[]);
  
  K = [];
  for p = 1 : num_te
      X1 = repmat(X(p,com(:,1)),num_tr,1);
      X1 = reshape(X1,1,[]);
      X2 = repmat(X(p,com(:,2)),num_tr,1);
      X2 = reshape(X2,1,[]);
      bd1 = repmat(bws(com(:,1)),num_tr,1);
      bd1 = reshape(bd1,1,[]);
      bd2 = repmat(bws(com(:,2)),num_tr,1);
      bd2 = reshape(bd2,1,[]);
      tmp = ((X1-Y1).^2)./bd1 + ((X2-Y2).^2)./bd2;
      K = [K; exp(-0.5*tmp)];
  end

%   length = 1;
%   for D_1 = 1:dimension-1
%       for D_2 = D_1+1:dimension
%           for i = 1:num_tr
%               for p =1:num_te
%                   K1(p,length) = exp(-0.5*((Y(i,D_1)-X(p,D_1))^2/bws(D_1)+...
%                       (Y(i,D_2)-X(p,D_2))^2/bws(D_2)));
%               end
%               length = length +1;
%           end
%       end
%   end

end



