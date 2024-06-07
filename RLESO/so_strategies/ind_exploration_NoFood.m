function newX_dec = ind_exploration_NoFood(X_dec,fitness,ind_fitness,C2,lb,ub)
    vec_flag = [1,-1] ;        
    [N,dim] = size(X_dec) ;
    newX_dec = zeros(1,dim);

    for j=1:1:dim
        rand_leader_index = floor(N*rand()+1);
        X_rand = X_dec(rand_leader_index, :);
        flag_index = floor(2*rand()+1);
        Flag=vec_flag(flag_index);
        A=exp(-fitness(rand_leader_index)/(ind_fitness+eps));
        newX_dec(1,j)=X_rand(j)+Flag*C2*A*((ub-lb)*rand+lb);
    end

end

