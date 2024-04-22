function newX_dec= GaussianMutation(X_dec,LB,UB)
    newX_dec = X_dec .*(1 + (UB-LB) .* randn(1));
end

