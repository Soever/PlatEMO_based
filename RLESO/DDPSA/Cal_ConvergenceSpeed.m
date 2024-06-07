function [S_G] = Cal_ConvergenceSpeed(X_now,X_pre,S_pre,c)
    d= size(X_now,2) ;
    X_diff = X_now - X_pre;
    nonzero_counts = sum(X_diff ~= 0, 2);
    R_D = nonzero_counts./d ; 
    sigma = std(X_diff,1,2); 
    S_G = (1-c) *S_pre +c * sum(sigma.*R_D) ;
end

 