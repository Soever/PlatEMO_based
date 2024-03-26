function [X,fitness] = Tent_Chaos(fitness,N,X,dim,Problem)
    avgF = mean(fitness);
    for j = 1:N
        if(fitness(j) < avgF)
            %% Cauchy mutation
            TempX_dec = X(j).dec.*(1 + tan(pi*(rand-0.5)));%eq.(28)
        else
            %% Tent-chaos
            TentZ0 = rand;
            TentValue =mod(2*TentZ0,1)+rand/(N*dim);%eq.(25)
            newX_dec=min(X.decs)+(max(X.decs)-min(X.decs))*TentValue;
            TempX_dec = (X(j).dec+newX_dec)./2;
        end
        TempX_dec = Problem.CalDec(TempX_dec) ;
        TempX=Problem.Evaluation(TempX_dec);
        if(TempX.obj<fitness(j))
            fitness(j)= TempX.obj;
            X(j) = TempX;
        end
    end
end

