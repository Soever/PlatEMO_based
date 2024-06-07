function X_new = DynamicOBL(Xmin_t,Xmax_t,X_best)
    r1 =rand() ;r2 = rand() ;
    X0 = r1*(Xmin_t+Xmax_t)-X_best ;
    X_new = r2*X0 + (1-r2) *X_best;
    
end

