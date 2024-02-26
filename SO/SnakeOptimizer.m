classdef SnakeOptimizer < ALGORITHM
    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            N = Problem.N;  % Population size
            Nm=round(N/2);
            Nf=N-Nm;
            
            T = Problem.maxFE / Problem.N;  % 最大迭代次数
            dim = Problem.D;  % 决策变量维度
            lb = Problem.lower;  % 下界
            ub = Problem.upper;  % 上界

            vec_flag = [1, -1];
            Threshold = 0.25;
            Thresold2 = 0.6;
            C1 = 0.5;
            C2 = 0.05;
            C3 = 2;
            X = Problem.Initialization();
            fitness = X.objs ;
            [GYbest, gbest] = min(fitness);
            Xfood = X(gbest);
            Xm=X(1:Nm);
            Xf=X(Nm+1:N);
            fitness_m=fitness(1:Nm);
            fitness_f=fitness(Nm+1:N);
            [fitnessBest_m, gbest1] = min(fitness_m);
            Xbest_m = Xm(gbest1);
            [fitnessBest_f, gbest2] = min(fitness_f);
            Xbest_f = Xf(gbest2);
            t=0 ;
            while Algorithm.NotTerminated(X)
            %   
                t = t+1 ;
                Temp=exp(-((t)/T));  %eq.(4)
                Q=C1*exp(((t-T)/(T)));%eq.(5)
                Xnewm = Xm.decs ;
                Xnewf =Xf.decs  ;
       
                if Q>1        Q=1;    end
                if Q<Threshold
                    for i=1:Nm
                        for j=1:1:dim
                            rand_leader_index = floor(Nm*rand()+1);
                            X_randm = Xm(rand_leader_index);
                            flag_index = floor(2*rand()+1);
                            Flag=vec_flag(flag_index); 
                            Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));%eq.(7)
                            Xnewm(i,j)=X_randm.dec(j)+Flag*C2*Am*((ub(j)-lb(j)*rand+lb(j)));%eq.(6)
                        end
                    end
                    for i=1:Nf
                        for j=1:1:dim
                            rand_leader_index = floor(Nf*rand()+1);
                            X_randf = Xf(rand_leader_index);
                            flag_index = floor(2*rand()+1);
                            Flag=vec_flag(flag_index);
                            Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));%eq.(9)
                            Xnewf(i,j)=X_randf.dec(j)+Flag*C2*Af*((ub(j)-lb(j))*rand+lb(j));%eq.(8)
                        end
                    end
                else %Exploitation Phase (Food Exists)
                    if Temp>Thresold2  %hot
                        for i=1:Nm
                            flag_index = floor(2*rand()+1);
                            Flag=vec_flag(flag_index);
                            for j=1:1:dim
                                Xnewm(i,j)=Xfood.dec(j)+C3*Flag*Temp*rand*(Xfood.dec(j)-Xnewm(i,j));%eq.(10)
                            end
                        end
                        for i=1:Nf
                            flag_index = floor(2*rand()+1);
                            Flag=vec_flag(flag_index);
                            for j=1:1:dim
                                Xnewf(i,j)=Xfood.dec(j)+Flag*C3*Temp*rand*(Xfood.dec(j)-Xnewf(i,j));%eq.(10)
                            end
                        end
                    else %cold
                        if rand>0.6 %fight
                            for i=1:Nm
                                for j=1:1:dim
                                    FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));%eq.(13)
                                    Xnewm(i,j)=Xnewm(i,j) +C3*FM*rand*(Q*Xbest_f.dec(j)-Xnewm(i,j));%eq.(11)
                                    
                                end
                            end
                            for i=1:Nf
                                for j=1:1:dim
                                    FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));%eq.(14)
                                    Xnewf(i,j)=Xnewf(i,j)+C3*FF*rand*(Q*Xbest_m.dec(j)-Xnewf(i,j));%eq.(12)
                                end
                            end
                        else%mating
                            for i=1:Nm
                                for j=1:1:dim
                                    Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));%eq.(17)
                                    Xnewm(i,j)=Xnewm(i,j) +C3*rand*Mm*(Q*Xnewf(i,j)-Xnewm(i,j));%eq.(15
                                end
                            end
                            for i=1:Nf
                                for j=1:1:dim
                                    Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));%eq.(18)
                                    Xnewf(i,j)=Xnewf(i,j) +C3*rand*Mf*(Q*Xnewm(i,j)-Xnewf(i,j));%eq.(16)
                                end
                            end
                            flag_index = floor(2*rand()+1);
                            egg=vec_flag(flag_index);
                            if egg==1;
                                [GYworst, gworst] = max(fitness_m);
                                Xnewm(gworst,:)=lb+rand*(ub-lb);%eq.(19)
                                [GYworst, gworst] = max(fitness_f);
                                Xnewf(gworst,:)=lb+rand*(ub-lb);%eq.(20)
                            end
                        end
                    end
                end
                for j=1:Nm
                    Xnewm = Problem.CalDec(Xnewm) ;
                    Xnewm_ind = Problem.Evaluation(Xnewm(j,:));
                    y = Xnewm_ind.obj ;
                    if y<fitness_m(j)
                        fitness_m(j)=y;
                        Xm(j)= Xnewm_ind;
                    end
                end
                
                [Ybest1,gbest1] = min(fitness_m);
                
                for j=1:Nf

                    Xnewf = Problem.CalDec(Xnewf) ;
                    Xnewf_ind =Problem.Evaluation(Xnewf(j,:)) ;
                    y = Xnewf_ind.obj ;
                    if y<fitness_f(j)
                        fitness_f(j)=y;
                        Xf(j)= Xnewf_ind;
                    end
                end
                X(1:Nm) = Xf ;
                X(Nm+1:end) = Xf ;
                [Ybest2,gbest2] = min(fitness_f);
                
                if Ybest1<fitnessBest_m
                    Xbest_m = Xm(gbest1);
                    fitnessBest_m=Ybest1;
                end
                if Ybest2<fitnessBest_f
                    Xbest_f = Xf(gbest2);
                    fitnessBest_f=Ybest2;
                    
                end
                if Ybest1<Ybest2
                    gbest_t(t)=min(Ybest1);
                else
                    gbest_t(t)=min(Ybest2);
                end
                if fitnessBest_m<fitnessBest_f
                    GYbest=fitnessBest_m;
                    Xfood=Xbest_m;
                else
                    GYbest=fitnessBest_f;
                    Xfood=Xbest_f;
                end

            end
        end
    end
end