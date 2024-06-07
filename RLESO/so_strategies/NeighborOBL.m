function Temp = NeighborOBL(k,Xbest,X_neighbor)
    ub = max(X_neighbor,[],1);
    lb = min(X_neighbor,[],1);
    Temp = (ub + lb)./2 + (ub + lb)./(2*k) - Xbest./k; %eq.(20)
    Flag4ub=Temp>ub;
    Flag4lb=Temp<lb;
    Temp=(Temp.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
end

