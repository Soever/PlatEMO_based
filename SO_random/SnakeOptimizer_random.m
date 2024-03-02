classdef SnakeOptimizer_random < ALGORITHM
    methods
        function main(Algorithm,Problem)
            global result_t ;
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
         
            result_t(1)=GYbest;
            X_personalbest = X ;
            t2_count =0 ;
            while Algorithm.NotTerminated(X)
            %   
                t = t+1 ;
                Temp=exp(-((t)/T));  %eq.(4)
                Q=C1*exp(((t-T)/(T)));%eq.(5)
                Xnewm = Xm.decs ;
                Xnewf =Xf.decs  ;
                if Q>1        Q=1;    end
                if Q<Threshold
                    %雄性探索
                    Xnewm = so_exploration(Nm,dim,Xm,C2,ub,lb);
                    Xnewf = so_exploration(Nf,dim,Xf,C2,ub,lb);
                else %Exploitation Phase (Food Exists)
                    if t2_count ==0 
                        X = X_personalbest ;
                        Xm=X(1:Nm);
                        Xf=X(Nm+1:N);
                        t2_count = t2_count+1 ;
                    end
                    if Temp>Thresold2  %hot
                        Xnewm =so_exploitation_food_hot(Nm,dim,Xfood,C3,Temp,Xm);
                        Xnewf =so_exploitation_food_hot(Nf,dim,Xfood,C3,Temp,Xf);
                    else %cold
                       [Xnewm,Xnewf] = so_exploitation_food_cold(Nm,Nf,dim,fitnessBest_m,fitnessBest_f,fitness_m,fitness_f,C3,Q,Xbest_m,Xbest_f,Xm,Xf,lb,ub);
                    end
                end
                if Q<Threshold
                    Xnewm = Problem.CalDec(Xnewm) ;
                    Xnewm = Problem.Evaluation(Xnewm);
                    fitness_m = Xnewm.objs ;
                    Xnewf = Problem.CalDec(Xnewf) ;
                    Xnewf = Problem.Evaluation(Xnewf);
                    fitness_f = Xnewf.objs ;
                    
                    Xm = Xnewm ;
                    Xf = Xnewf ;
                    X(1:Nm) = Xm ;
                    X(Nm+1:N) = Xf ;
                    for ii = 1:N
                        if X(ii).obj < X_personalbest(ii).obj
                            X_personalbest(ii) = X(ii) ;
                        end
                    end

                else
                    %更新雄性种群
                    for j=1:Nm
                        Xnewm = Problem.CalDec(Xnewm) ;
                        Xnewm_ind = Problem.Evaluation(Xnewm(j,:));
                        y = Xnewm_ind.obj ;
                        if y<fitness_m(j)
                            fitness_m(j)=y;
                            Xm(j)= Xnewm_ind;
                        end
                    end
                    %更新雌性种群
                    for j=1:Nf
                        Xnewf = Problem.CalDec(Xnewf) ;
                        Xnewf_ind =Problem.Evaluation(Xnewf(j,:)) ;
                        y = Xnewf_ind.obj ;
                        if y<fitness_f(j)
                            fitness_f(j)=y;
                            Xf(j)= Xnewf_ind;
                        end
                    end 
                end    
                X(1:Nm) = Xm ;
                X(Nm+1:end) = Xf ;
                [Ybest1,gbest1] = min(fitness_m);
                [Ybest2,gbest2] = min(fitness_f);
                %更新雄性、雌性最优
                if Ybest1<fitnessBest_m
                    Xbest_m = Xm(gbest1);
                    fitnessBest_m=Ybest1;
                end
                if Ybest2<fitnessBest_f
                    Xbest_f = Xf(gbest2);
                    fitnessBest_f=Ybest2;
                end
                %记录全局最优
                if Ybest1<Ybest2
                    gbest_t(t)=min(Ybest1);
                else
                    gbest_t(t)=min(Ybest2);
                end
                if t == T-1
                    result_t(2:t+1) = gbest_t ;
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