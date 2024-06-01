function newX_dec = FastRandomOBL(X_dec,LB,UB)
    dim = size(X_dec,2);
    newX_dec = zeros(1,dim) ;
    M = (LB+UB)/2 ;
    for i = 1:dim
        if abs(X_dec(i)) < abs(M)
            newX_dec(1,i) = M(1,i)+rand^2*sin(2*pi*rand)*X_dec(1,i)/2 ;
        else
            newX_dec(1,i) = M(1,i)-rand^2*sin(2*pi*rand)*X_dec(1,i)/2 ;
        end
    end
end