function newX_dec= GaussianMutation(X_dec,LB,UB)
[N,dim] = size(X_dec);
    newX_dec = zeros(N,dim) ;
    for i=1:N
        newX_dec(i,:) = X_dec(i,:) .*(1 + (UB-LB) .* randn(1));
    end    
end

