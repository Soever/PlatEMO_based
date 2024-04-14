function newX_dec = exploration_NoFood(X_dec,fitness,C2,lb,ub)
    vec_flag = [1,-1] ;        
    [N,dim] = size(X_dec) ;
    newX_dec = zeros(N,dim);
    for i=1:N
        for j=1:1:dim
            rand_leader_index = floor(N*rand()+1);
            X_rand = X_dec(rand_leader_index, :);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            A=exp(-fitness(rand_leader_index)/(fitness(i)+eps));
            newX_dec(i,j)=X_rand(j)+Flag*C2*A*((ub-lb)*rand+lb);
        end
    end
end

