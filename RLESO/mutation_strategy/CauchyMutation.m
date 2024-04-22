function ind_dec = CauchyMutation(ind_dec)
    [N,dim] = size(ind_dec);
    newX_dec = zeros(N,dim) ;
    for i=1:N
        newX_dec(i,:) = ind_dec(i,:) .* (1 + tan(pi*(rand-0.5))) ;
    end

end

