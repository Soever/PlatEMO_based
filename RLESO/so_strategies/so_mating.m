function [newXm_dec, newXf_dec]= so_mating(Xm_dec,Xf_dec,fitness_m,fitness_f,C3,Q,lb,ub)
       
    [Nm,dim] = size(Xm_dec) ;
    newXm_dec = zeros(Nm,dim);    
    for i=1:Nm
        for j=1:1:dim
            Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));
            newXm_dec(i,j)=Xm_dec(i,j) +C3*rand*Mm*(Q*Xf_dec(i,j)-Xm_dec(i,j));%eq.(10)
        end
    end
    [Nf,dim] = size(Xf_dec) ;
    newXf_dec = zeros(Nf,dim);  
    for i=1:Nf
        for j=1:1:dim
            Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));
            newXf_dec(i,j)=Xf_dec(i,j) +C3*rand*Mf*(Q*Xm_dec(i,j)-Xf_dec(i,j));%eq.(11)
        end
    end
    % [~, index]=sort(fitness_m);
    % [~, index1]= sort(fitness_f);%排序
    % newXm_dec(index(end-3),:)=lb+rand*(ub-lb);
    % newXm_dec(index(end-2),:)=lb+rand*(ub-lb);
    % newXm_dec(index(end-1),:)=lb+rand*(ub-lb);
    % newXm_dec(index(end),:)=lb+rand*(ub-lb);
    % newXf_dec(index1(end-3),:)=lb+rand*(ub-lb);
    % newXf_dec(index1(end-1),:)=lb+rand*(ub-lb);
    % newXf_dec(index1(end-2),:)=lb+rand*(ub-lb);
    % newXf_dec(index1(end),:)=lb+rand*(ub-lb);
end

