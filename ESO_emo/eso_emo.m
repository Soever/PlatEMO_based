classdef eso_emo < ALGORITHM
    methods
        function main(Algorithm,Problem)

        %[Threshold,Thresold2,C1,C2,C3] = Algorithm.ParameterSet(0.25,0.6,0.5,0.05,2);
        % [Xfood, fval,gbest_t,Trajectories,fitness_history, position_history] = 
        % ESO(N,T,lb,ub,dim,fobj)
        % N 种群个数 T 总迭代次数 lb 
        %Trajectories=zeros(N,T);
        %% initial
        
        Threshold=0.25;
        Thresold2= 0.6;
        T = Problem.maxFE / Problem.N; 
        N = Problem.N;
        dim =  Problem.D;
        C1=0.5*ones(1,T);
        C2=0.05*ones(1,T);
        C3=2*ones(1,T);
        %%
        [C1,C2,C3] = DynamicUpdate(C1,C2,C3) ;
        %%
        % for l= 1:T
        %     C1(1,l)=C1(1,l)+0.1*sin((pi/2)*((rand).^4));%eq.(14)
        %     C2(1,l)=C2(1,l)+0.001*cos((pi/2)*((rand).^4));%eq.(15)
        %     C3(1,l)=C3(1,l)-2*sin(0.5*pi*(l/T)^4);%eq.(16)
        % end



        lb = Problem.lower; 
        ub = Problem.upper;

        
        %% X=initializationNew(N,dim,ub,lb,fobj);
        X = Problem.Initialization();


        
        t1=ones(1,T);
        t2=ones(1,T);
        a=0.05;
        %%
        [t1,t2] = SineCosineCompositePerturbationFactors(t1,t2,a);
        %%
        % for l= 1:T
        %     t1(1,l)=1+(0.0001*(sin(a*4*pi*l)+cos(a*6*pi*l)))*exp((pi/100)*(0.25*(T-l)));
        %     t2(1,l)=1+(0.0001*(cos(a*4*pi*l)+sin(a*6*pi*l)))*exp((pi/100)*(0.25*(T-l)));
        % end
        [GYbest, gbest] = min(X.objs);
        Xfood = X(gbest);
        %% Diving the swarm into two equal groups males and females
        Nm=round(N/2);%eq.(2&3)
        Nf=N-Nm;
        Xm=X(1:Nm);
        Xf=X(Nm+1:N);
 
        [fitnessBest_m, gbest1] = min(Xm.objs);Xbest_m = Xm(gbest1);
        [fitnessBest_f, gbest2] = min(Xf.objs);Xbest_f = Xf(gbest2);
        t=0 ;
        %% Main loop
        while Algorithm.NotTerminated(X)
            t = t+1;
             %% the principle of convex lens imaging
            % k = 10*(1-2*(t/T)^2);% scaling factor eq.(19)
            % 
            % TempXm = convex_lens_imaging(ub,lb,k,Xbest_m.dec,Problem);
            % if(TempXm.obj<GYbest)
            %     fitnessBest_m=TempXm.obj ;
            %     Xbest_m = TempXm;
            %     Xm(gbest1) = TempXm;
            % end
            % 
            % TempXf = convex_lens_imaging(ub,lb,k,Xbest_f.dec,Problem);
            % if(TempXf.obj<GYbest)
            %     fitnessBest_f=TempXf.obj ;
            %     Xbest_f = TempXf;
            %     Xf(gbest2) = TempXf;
            % end
            k = 10*(1-2*(t/T)^2);% scaling factor eq.(19)
            TempXm_dec = ConvexLensImaging(k,Xbest_m.dec,ub,ub1,lb,lb1);
            TempXm=Problem.Evaluation(TempXm_dec);
          
            if(TempXm.obj<GYbest)
                fitnessBest_m=TempXm.obj ;
                Xbest_m = TempXm;
                Xm(gbest1) = TempXm;
            end
            TempXf_dec = ConvexLensImaging(k,Xbest_f.dec,ub,ub1,lb,lb1);
            TempXf=Problem.Evaluation(TempXf_dec);
            if(fitTemp<GYbest)
                fitnessBest_f=TempXf.obj ;
                Xbest_f = TempXf;
                Xf(gbest2) = TempXf;
            end
   
            Temp=exp(-((t)/T));  %eq.(4)
            Q=C1(1,t)*exp(((t-T)/(T)));%eq.(5)
            if Q>1        Q=1;    end
            %% Exploration Phase (No Food)
            if Q<Threshold
                %探索
                Xnewm_dec = so_exploration(Nm,dim,Xm,C2(1,t),ub,lb);
                Xnewf_dec = so_exploration(Nf,dim,Xf,C2(1,t),ub,lb);
              
            %% Exploitation Phase (Food Exists)
            else
                if Temp>Thresold2  %hot
                    Xnewm_dec =so_exploitation_food_hot(Nm,dim,Xfood,C3(1,t),Temp,Xm);
                    Xnewf_dec =so_exploitation_food_hot(Nf,dim,Xfood,C3(1,t),Temp,Xf);
                else %cold
                   [Xnewm_dec,Xnewf_dec] = improved_exploitation_food_cold(t1(1,t),t2(1,t),Nm,Nf,dim,C3(1,t),Q,Xbest_m,Xbest_f,Xm,Xf,lb,ub);
                end
            end
        
            %% 更新雄性种群
            Xm = CalNewPop(Xm,Nm,Xnewm_dec,Problem);
            Xf = CalNewPop(Xf,Nf,Xnewf_dec,Problem);
            %% 混沌与柯西
            [Xm,fitness_m] = Tent_Chaos(Nm,Xm,dim,Problem);
            [Xf,fitness_f] = Tent_Chaos(Nf,Xf,dim,Problem);
        
            [Ybest1,gbest1] = min(Xm.objs);
            [Ybest2,gbest2] = min(Xf.objs);
            X(1:Nm) = Xm ;
            X(Nm+1:end) = Xf ;
            if Xm(gbest1).obj<Xbest_m.obj
                Xbest_m = Xm(gbest1);
            end
            if Xf(gbest2).obj<Xbest_f.obj
                Xbest_f = Xf(gbest2);
            end
            if Xbest_m.obj < Xbest_f.obj 
                Xfood=Xbest_m;
            else
                Xfood=Xbest_f;
            end
            gbest_t(t)=min(Ybest1,Ybest2);%第t代的最优值
            


        end
            fval = Xfood.obj;
        end
    end
    
end

function X =CalNewPop(X,N,Xnew_dec,Problem)
    fitness=X.objs;
    Xnew_dec = Problem.CalDec(Xnew_dec) ;
    Xnew = Problem.Evaluation(Xnew_dec);
    for j=1:N
        if Xnew(j).obj<fitness(j)
            fitness(j)=Xnew(j).obj;
            X(j)= Xnew(j);
        end
    end

end


