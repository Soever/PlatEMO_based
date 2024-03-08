classdef SnakeOptimizer_pro < ALGORITHM
    methods
        function main(Algorithm,Problem)
            global result_t ;
            %% Parameter setting
            N = Problem.N;  % Population size
            Nm=round(N/2);
            Nf=N-Nm;
            
            T = Problem.maxFE / (Problem.N*2);  % 最大迭代次数
            dim = Problem.D;  % 决策变量维度
            lb = Problem.lower;  % 下界
            ub = Problem.upper;  % 上界

            vec_flag = [1, -1];
            Threshold = 0.25;
            Thresold2 = 0.6;
            C1 = 0.5;
            C2 = 0.05;
            C3 = 2;
            RR = 0.7 ;
            X = Problem.Initialization();
            
            for i = 1:N
                Schedule(i,1:dim) = unifrnd(0.9,1,1,dim);
            end
            fitness = X.objs ;
            [GYbest, gbest] = min(fitness);
            Xfood = X(gbest);
            Xm=X(1:Nm);
            Xf=X(Nm+1:N);
            [fitnessBest_m, gbest1] = min(Xm.objs);
            Xbest_m = Xm(gbest1);
            [fitnessBest_f, gbest2] = min(Xf.objs);
            Xbest_f = Xf(gbest2);
            t=0 ;
            result_t(1)=GYbest;
            while Algorithm.NotTerminated(X)
            %   
                t = t+1 ;
                Temp=exp(-((t)/T));  %eq.(4)
                Q=C1*exp(((t-T)/(T)));%eq.(5)
                if Q>1        Q=1;    end
                if Q<Threshold
                    %雄性探索
                    Xnewm = so_exploration(Nm,dim,Xm,C2,ub,lb);
                    Xnewf = so_exploration(Nf,dim,Xf,C2,ub,lb);

                else %Exploitation Phase (Food Exists)
                    if Temp>Thresold2  %hot
                        Xnewm =so_exploitation_food_hot(Nm,dim,Xfood,C3,Temp,Xm);
                        Xnewf =so_exploitation_food_hot(Nf,dim,Xfood,C3,Temp,Xf);
                    else %cold
                       [Xnewm,Xnewf] = so_exploitation_food_cold(Nm,Nf,dim,C3,Q,Xbest_m,Xbest_f,Xm,Xf,lb,ub);
                    end
                end
 
                

                %更新雄性种群
                for j=1:Nm
                    Xnewm = Problem.CalDec(Xnewm) ;
                    Xnewm_ind = Problem.Evaluation(Xnewm(j,:));
                    if Xnewm_ind.obj<Xm(j).obj 
                        Xm(j)= Xnewm_ind ;
                        Schedule(j) = Schedule(j) + Schedule(j) * (RR);
                    else
                        Schedule(j) = Schedule(j) - Schedule(j) * (RR);
                    end
                end
                
                %更新雌性种群
                for j=1:Nf
                    Xnewf = Problem.CalDec(Xnewf) ;
                    Xnewf_ind =Problem.Evaluation(Xnewf(j,:)) ;
                    if Xnewf_ind.obj<Xf(j).obj
                        Xf(j)= Xnewf_ind;
                        Schedule(j+Nm) = Schedule(j+Nm) + Schedule(j+Nm) * (RR);
                    else
                        Schedule(j+Nm) = Schedule(j+Nm) - Schedule(j+Nm) * (RR);
                    end
                end
                [Ybest1,gbest1] = min(Xm.objs);
                [Ybest2,gbest2] = min(Xf.objs);
                X(1:Nm) = Xm ;
                X(Nm+1:end) = Xf ;
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

                Tau = (t/T);
                Xnew = pro_method(N,dim,X,Xfood,Tau,Schedule) ;
                Xnew = Problem.CalDec(Xnew) ;
                Xnew = Problem.Evaluation(Xnew);
                for i = 1:N
                    if Xnew(i).obj < X(i).obj
                        X(i) = Xnew(i);
                        Schedule(j) = Schedule(j) + Schedule(j) * (RR);
                    else
                        Schedule(j) = Schedule(j) - Schedule(j) * (RR);
                    end
                end
                Xm=X(1:Nm);
                Xf=X(Nm+1:N);
                [Ybest1,gbest1] = min(Xm.objs);
                [Ybest2,gbest2] = min(Xf.objs);
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