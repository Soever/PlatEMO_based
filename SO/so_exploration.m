function Xnew = so_exploration(N,dim,X,fitness,C2,Xnew,ub,lb)
vec_flag = [1, -1];
    for i=1:N
        for j=1:1:dim
            rand_leader_index = floor(N*rand()+1);
            X_randm = X(rand_leader_index);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index); 
            Am=exp(-fitness(rand_leader_index)/(fitness(i)+eps));%eq.(7)
            Xnew(i,j)=X_randm.dec(j)+Flag*C2*Am*((ub(j)-lb(j)*rand+lb(j)));%eq.(6)
        end
    end
end