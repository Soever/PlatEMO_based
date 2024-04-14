
function [Xfood, fval,gbest_t,Trajectories,fitness_history, position_history] = RLSO_FitnessDiversityGrid(N,T,lb,ub,dim,fobj)


%% initial

Threshold=0.25;
Thresold2= 0.6;
C1=0.5*ones(1,T);
C2=0.05*ones(1,T);
C3=2*ones(1,T);
stage = [1,2,3,4] ;
t1=ones(1,T);
t2=ones(1,T);
ub1 = ub.*ones(1,dim);
lb1 = lb.*ones(1,dim);
if(max(size(ub)) == 1)
    ub2 = ub.*ones(1,dim);
    lb2 = lb.*ones(1,dim);
end

%% X=initializationNew(N,dim,ub,lb,fobj);
X=lb+rand(N,dim)*(ub-lb);%eq.(1)
fitness=zeros(1,N);
for i=1:N
    fitness(i)=feval(fobj,X(i,:));
end
Trajectories=zeros(N,T);
position_history=zeros(N,T,dim);
fitness_history=zeros(N,T);

gbest_t = zeros(1,T);
[GYbest, gbest] = min(fitness);
Xfood = X(gbest,:);
%% Diving the swarm into two equal groups males and females
Nm=round(N/2);%eq.(2&3)
Nf=N-Nm;
Xm=X(1:Nm,:);
Xf=X(Nm+1:N,:);
fitness_m=fitness(1:Nm);fitness_f=fitness(Nm+1:N);
[fitnessBest_m, gbest1] = min(fitness_m);Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);Xbest_f = Xf(gbest2,:);

%% Main loop
for t = 1:T
    Positions=[Xm;Xf];
    for i=1:size(Positions,1)
        position_history(i,t,:)=Positions(i,:);
        Trajectories(:,t)=Positions(:,1);
        fitness_history(i,t)=fobj(Positions(i,:));
    end
    Temp=exp(-((t)/T));  %eq.(4)
    Q=C1(1,t)*exp(((t-T)/(T)));%eq.(5)
    if Q>1        Q=1;    end
    if Q<Threshold
         newXm_dec = exploration_NoFood(Xm,fitness_m,C2(1,t),lb,ub);
         newXf_dec = exploration_NoFood(Xf,fitness_f,C2(1,t),lb,ub);
    else
        if Temp>Thresold2
            newXm_dec = exploit_Food(Xm,Xfood,Temp,C3(1,t));
            newXf_dec = exploit_Food(Xf,Xfood,Temp,C3(1,t));
        else 
            if rand>0.6 
                newXm_dec = so_fight(Xm,fitness_m,Xbest_f,fitnessBest_f,t1(1,t),C3(1,t),Q) ;
                newXf_dec = so_fight(Xf,fitness_f,Xbest_m,fitnessBest_m,t1(1,t),C3(1,t),Q) ;
            else
                [newXm_dec, newXf_dec] = so_mating(Xm,Xf,fitness_m,fitness_f,C3(1,t),Q,lb,ub);
            end
        end
    end
    [Xm,fitness_m] = Evaluation(Xm,newXm_dec,fitness_m,lb,ub,fobj);
    [Xf,fitness_f] = Evaluation(Xf,newXf_dec,fitness_f,lb,ub,fobj);
    [Xbest_m,Xbest_f,fitness_m,fitness_f,GYbest,Xfood] = updateXbest(Xm,Xf,fitness_m,fitness_f,Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f);
    gbest_t(1,t) = GYbest ;
end
    fval = GYbest;
end

function q = updataQtable(s,a,r,q)
    s_next = a;
    [q_target,index] = max(q(s_next,:));
    q(s,a) =q(s,a)+0.1*(r+0.9*q_target-q(s,a)) ;
    
end





function [X,fitness] = Evaluation(X,newX_dec,fitness,lb,ub,fobj)
    N = size(X,1);
    for j=1:N
        Flag4ub=newX_dec(j,:)>ub;
        Flag4lb=newX_dec(j,:)<lb;
        newX_dec(j,:)=(newX_dec(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,newX_dec(j,:));
        if y<fitness(j)
            fitness(j)=y;
            X(j,:)= newX_dec(j,:);
        end
    end
end

function [Xbest_m,Xbest_f,fitness_m,fitness_f,GYbest,Xfood] = updateXbest(Xm,Xf,fitness_m,fitness_f,Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f)
    [Ybest1,gbest1] = min(fitness_m);
    [Ybest2,gbest2] = min(fitness_f);
    if Ybest1<fitnessBest_m
        Xbest_m = Xm(gbest1,:);
        fitnessBest_m=Ybest1;
    end
    if Ybest2<fitnessBest_f
        Xbest_f = Xf(gbest2,:);
        fitnessBest_f=Ybest2;
    end
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        Xfood=Xbest_m;
    else
        GYbest=fitnessBest_f;
        Xfood=Xbest_f;
    end
end