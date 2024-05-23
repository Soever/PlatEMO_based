function newX_dec= best_worst_position(X_dec,Xbest_dec,Xworst_dec,t1,C3,Q)
       
    [N,dim] = size(X_dec) ;
    newX_dec = zeros(N,dim);    
    for i=1:N
        for j=1:1:dim
            newX_dec(i,j)=t1*X_dec(i,j) +rand*(Q*Xbest_dec(j)-X_dec(i,j)) -rand*(Q*Xworst_dec(j)-X_dec(i,j));%eq.(8)
        end
    end
end

