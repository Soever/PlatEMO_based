classdef ESO2 < ALGORITHM
    methods
        function main(Algorithm,Problem)

        %[Threshold,Thresold2,C1,C2,C3] = Algorithm.ParameterSet(0.25,0.6,0.5,0.05,2);
        % [Xfood, fval,gbest_t,Trajectories,fitness_history, position_history] = 
        % ESO(N,T,lb,ub,dim,fobj)
        % N 种群个数 T 总迭代次数 lb 
        %Trajectories=zeros(N,T);
        %% initial
        vec_flag=[1,-1];
        Threshold=0.25;
        Thresold2= 0.6;
        T = Problem.maxFE / Problem.N; 
        N = Problem.N;
        dim =  Problem.D;
        C3=zeros(1,T);
        for l= 1:T
            C3(1,l)=2;%eq.(16)
        end


        lb = Problem.lower; 
        ub = Problem.upper;

        
        %% X=initializationNew(N,dim,ub,lb,fobj);
        X = Problem.Initialization();
        fitness = X.objs ;

        
        t1=zeros(1,T);
        t2=zeros(1,T);
        a=0.05;
        for l= 1:T
            t1(1,l)=1;
            t2(1,l)=1;
        end
        [GYbest, gbest] = min(fitness);
        Xfood = X(gbest);
        %% Diving the swarm into two equal groups males and females
        Nm=round(N/2);%eq.(2&3)
        Nf=N-Nm;
        Xm=X(1:Nm);
        Xf=X(Nm+1:N);
        fitness_m=fitness(1:Nm);
        fitness_f=fitness(Nm+1:N);
        [fitnessBest_m, gbest1] = min(fitness_m);Xbest_m = Xm(gbest1);
        [fitnessBest_f, gbest2] = min(fitness_f);Xbest_f = Xf(gbest2);
        t=0 ;
        %% Main loop
        while Algorithm.NotTerminated(X)
            t = t+1;
             %% the principle of convex lens imaging
            k = 10*(1-2*(t/T)^2);% scaling factor eq.(19)
            
            TempXm = convex_lens_imaging(ub,lb,k,Xbest_m.dec,Problem);
            if(TempXm.obj<GYbest)
                fitnessBest_m=TempXm.obj ;
                Xbest_m = TempXm;
                Xm(gbest1) = TempXm;
            end

            TempXf = convex_lens_imaging(ub,lb,k,Xbest_f.dec,Problem);
            if(TempXf.obj<GYbest)
                fitnessBest_f=TempXf.obj ;
                Xbest_f = TempXf;
                Xf(gbest2) = TempXf;
            end
        
            C1=0.5;%eq.(14)
            C2=0.05;%eq.(15)
            Temp=exp(-((t)/T));  %eq.(4)
            Q=C1*exp(((t-T)/(T)));%eq.(5)
            if Q>1        Q=1;    end
            %% Exploration Phase (No Food)
            if Q<Threshold
                %探索
                Xnewm = so_exploration(Nm,dim,Xm,C2,ub,lb);
                Xnewf = so_exploration(Nf,dim,Xf,C2,ub,lb);
              
            %% Exploitation Phase (Food Exists)
            else
                if Temp>Thresold2  %hot
                    if Temp>Thresold2  %hot
                        Xnewm =so_exploitation_food_hot(Nm,dim,Xfood,C3(1,t),Temp,Xm);
                        Xnewf =so_exploitation_food_hot(Nf,dim,Xfood,C3(1,t),Temp,Xf);
                    else %cold
                       [Xnewm,Xnewf] = improved_exploitation_food_cold(t1(1,t),t2(1,t),Nm,Nf,dim,C3(1,t),Q,Xbest_m,Xbest_f,Xm,Xf,lb,ub);
                    end
                end
            end
        
            %% 更新雄性种群
            [Xm,fitness_m] = CalNewPop(Xm,Nm,Xnewm,fitness_m,Problem);
            [Xf,fitness_f] = CalNewPop(Xf,Nf,Xnewf,fitness_f,Problem);
            %% 混沌与柯西
            % [Xm,fitness_m] = Tent_Chaos(fitness_m,Nm,Xm,dim,Problem);
            % [Xf,fitness_f] = Tent_Chaos(fitness_f,Nf,Xf,dim,Problem);
        
            [Ybest1,gbest1] = min(fitness_m);
            [Ybest2,gbest2] = min(fitness_f);
            if Ybest1<fitnessBest_m
                Xbest_m = Xm(gbest1);
                fitnessBest_m=Ybest1;
            end
            if Ybest2<fitnessBest_f
                Xbest_f = Xf(gbest2);
                fitnessBest_f=Ybest2;
            end

            gbest_t(t)=min(Ybest1,Ybest2);%第t代的最优值
            

            if fitnessBest_m<fitnessBest_f
                GYbest=fitnessBest_m;
                Xfood=Xbest_m;
            else
                GYbest=fitnessBest_f;
                Xfood=Xbest_f;
            end
        end
            fval = GYbest;
        end
    end
    
end

function [X,fitness] =CalNewPop(X,N,Xnew_dec,fitness,Problem)
    Xnew_dec = Problem.CalDec(Xnew_dec) ;
    Xnew = Problem.Evaluation(Xnew_dec);
    for j=1:N
        if Xnew(j).obj<fitness(j)
            fitness(j)=Xnew(j).obj;
            X(j)= Xnew(j);
        end
    end
end


