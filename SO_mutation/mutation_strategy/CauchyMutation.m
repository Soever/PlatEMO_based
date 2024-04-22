function ind_dec = CauchyMutation(ind_dec)
    ind_dec = ind_dec .* (1 + tan(pi*(rand-0.5))) ;
end

