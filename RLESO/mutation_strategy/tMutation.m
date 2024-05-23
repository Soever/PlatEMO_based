function newX_dec= tMutation(X_dec,df)
    [N,dim] = size(X_dec);
    newX_dec = zeros(N,dim) ;
    for i=1:N
        t_randoms = trnd(df, 1);
        newX_dec(i,:) = X_dec(i,:) .* (1 + t_randoms);
    end   
end

