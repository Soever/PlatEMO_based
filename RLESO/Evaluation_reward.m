function [newX,fitness,reward] = Evaluation_reward(X,newX_dec,fitness,lb,ub,fobj)
    N = size(X,1);
    fitness_new = zeros(size(fitness));
    fitness_old = fitness ;
    for j=1:N
        Flag4ub=newX_dec(j,:)>ub;
        Flag4lb=newX_dec(j,:)<lb;
        newX_dec(j,:)=(newX_dec(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness_new(j) = feval(fobj,newX_dec(j,:));
        % 择优
        if fitness_new(j)  < fitness(j)
            X(j,:) = newX_dec(j,:) ;
            fitness(j) = fitness_new(j);
        end
    end
    
    if min(fitness_new) <min(fitness_old) 
        if mean(fitness_new) <mean(fitness_old)
            reward = 3;
        else
            reward =1 ;
        end
    else
        if mean(fitness_new) <mean(fitness_old)
            reward = -1;
        else
            reward =-3 ;
        end
    end

    % 不择优
    newX = X;
    % fitness = fitness_new;
end