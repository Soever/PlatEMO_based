function newX_dec= so_fight(X_dec,fitness,Xbest_dec,fitnessBest,t1,C3,Q)
       
    [N,dim] = size(X_dec) ;
    newX_dec = zeros(N,dim);    
    for i=1:N
        for j=1:1:dim
            F=exp(-(fitnessBest)/(fitness(i)+eps));
            newX_dec(i,j)=t1*X_dec(i,j) +C3*F*rand*(Q*Xbest_dec(j)-X_dec(i,j));%eq.(8)
        end
    end
end

