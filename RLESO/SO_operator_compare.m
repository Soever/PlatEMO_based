function [Xfood, fval,gbest_t,Trajectories,fitness_history, position_history,history] = SO_operator_compare(N,T,lb,ub,dim,fobj,stage)
%initial 
vec_flag=[1,-1];
Threshold=0.25;
Thresold2= 0.6;
% C1=0.5;
% C2=.05;
% C3=2;
C1=0.5*ones(1,T);
C2=0.05*ones(1,T);
C3=2*ones(1,T);
t1=ones(1,T);
t2=ones(1,T);
history = zeros(1,T);
X=lb+rand(N,dim)*(ub-lb);
for i=1:N
 fitness(i)=feval(fobj,X(i,:));   
end
Trajectories=zeros(N,T);
position_history=zeros(N,T,dim);
fitness_history=zeros(N,T);
[GYbest, gbest] = min(fitness);
Xfood = X(gbest,:);
%Diving the swarm into two equal groups males and females
Nm=round(N/2);%eq.(2&3)
Nf=N-Nm;
Xm=X(1:Nm,:);
Xf=X(Nm+1:N,:);
fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:N);
[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
for t = 1:T
  Temp=exp(-((t)/T));  %eq.(4)
  Q=C1(1,t)*exp(((t-T)/(T)));%eq.(5)
  Positions=[Xm;Xf];
    for i=1:size(Positions,1)
        position_history(i,t,:)=Positions(i,:);
        Trajectories(:,t)=Positions(:,1);
        fitness_history(i,t)=fobj(Positions(i,:));
    end
    if Q>1        Q=1;    end
    if stage == 1
        newXm_dec = exploration_NoFood(Xm,fitness_m,C2(1,t),lb,ub);
        newXf_dec = exploration_NoFood(Xf,fitness_f,C2(1,t),lb,ub);
    elseif stage == 2
        newXm_dec = exploit_Food(Xm,Xfood,Temp,C3(1,t));
        newXf_dec = exploit_Food(Xf,Xfood,Temp,C3(1,t));
    elseif stage == 3
        newXm_dec = so_fight(Xm,fitness_m,Xbest_f,fitnessBest_f,t1(1,t),C3(1,t),Q) ;
        newXf_dec = so_fight(Xf,fitness_f,Xbest_m,fitnessBest_m,t1(1,t),C3(1,t),Q) ;
    else
        [newXm_dec,newXf_dec] = so_mating(Xm,Xf,fitness_m,fitness_f,C3(1,t),Q,lb,ub);
    end
    % if Temp>Thresold2  %hot
        % newXm_dec = exploit_Food(Xm,Xfood,Temp,C3(1,t));
        % newXf_dec = exploit_Food(Xf,Xfood,Temp,C3(1,t));
    % else %cold
    %     if rand>0.6
    %         newXm_dec = so_fight(Xm,fitness_m,Xbest_f,fitnessBest_f,t1(1,t),C3(1,t),Q) ;
    %         newXf_dec = so_fight(Xf,fitness_f,Xbest_m,fitnessBest_m,t1(1,t),C3(1,t),Q) ;
    %     else
    %         [newXm_dec,newXf_dec] = so_mating(Xm,Xf,fitness_m,fitness_f,C3(1,t),Q,lb,ub);
    %     end
    % end
   
    Xnewm = newXm_dec ;
    Xnewf = newXf_dec ;
    hm = zeros(1,Nm);
    for j=1:Nm
         Flag4ub=Xnewm(j,:)>ub;
         Flag4lb=Xnewm(j,:)<lb;
        Xnewm(j,:)=(Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewm(j,:));
         hm(j) = y ;
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
        end
    end
    
    [Ybest1,gbest1] = min(fitness_m);
    hf = zeros(1,Nm);
    for j=1:Nf
         Flag4ub=Xnewf(j,:)>ub;
         Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewf(j,:));
        hf(j) = y ;
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);
        end
    end
    history(1,t) = min([min(hm),min(hf)]);
    [Ybest2,gbest2] = min(fitness_f);
    
    if Ybest1<fitnessBest_m
        Xbest_m = Xm(gbest1,:);
        fitnessBest_m=Ybest1;
    end
    if Ybest2<fitnessBest_f
        Xbest_f = Xf(gbest2,:);
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
fval = GYbest;
end





