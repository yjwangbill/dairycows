function K = espKernels_group(X, Y, bws, group)
% Computes Kernels using elementary symmetric polynomials.
% X :testMatrix, Y : training matrix
% return K num_te*(num_tr * (dimension-1)dimension/2 )

  % prelims
  [num_tr, dimension] = size(Y);
  num_te = size(X, 1);
  K = [];
  
  for g = 1 : length(group)
      k = group{g};
      for p = 1 : num_tr  
          Xk = X(:,k);
          Yk = repmat(Y(p,k),num_te,1);
          bdk = repmat(bws(k),num_te,1);
          tmpsum = sum(((Xk-Yk).^2)./bdk,2);
          K = [K exp(-0.5*tmpsum)];
      end
  end

end



