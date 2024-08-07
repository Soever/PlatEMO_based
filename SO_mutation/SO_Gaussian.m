
function [Xfood, fval,gbest_t,Trajectories,fitness_history, position_history] = SO_Gaussian(N,T,lb,ub,dim,fobj)


%% initial
vec_flag=[1,-1];
Threshold=0.25;
Thresold2= 0.6;
C1=0.5*ones(1,T);
C2=0.05*ones(1,T);
C3=2*ones(1,T);

t1=ones(1,T);
t2=ones(1,T);
a=0.05;


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


[GYbest, gbest] = min(fitness);
Xfood = X(gbest,:);
%% Diving the swarm into two equal groups males and females
Nm=round(N/2);%eq.(2&3)
Nf=N-Nm;
Xm=X(1:Nm,:);
Xf=X(Nm+1:N,:);
fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:N);
[fitnessBest_m, gbest1] = min(fitness_m);
 BestIndex1=gbest1;
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
 BestIndex2=gbest2;
Xbest_f = Xf(gbest2,:);

%% Main loop
for t = 1:T
        Positions=[Xm;Xf];
    for i=1:size(Positions,1)
        position_history(i,t,:)=Positions(i,:);
        Trajectories(:,t)=Positions(:,1);
        fitness_history(i,t)=fobj(Positions(i,:));
    end
     %% the principle of convex lens imaging
   

    Temp=exp(-((t)/T));  %eq.(4)
    Q=C1(1,t)*exp(((t-T)/(T)));%eq.(5)
    if Q>1        Q=1;    end
    %% Exploration Phase (No Food)
    if Q<Threshold
        for i=1:Nm
            for j=1:1:dim
                rand_leader_index = floor(Nm*rand()+1);
                X_randm = Xm(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));
                Xnewm(i,j)=X_randm(j)+Flag*C2(1,t)*Am*((ub-lb)*rand+lb);%eq.(5)
            end
        end
        for i=1:Nf
            for j=1:1:dim
                rand_leader_index = floor(Nf*rand()+1);
                X_randf = Xf(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));
                Xnewf(i,j)=X_randf(j)+Flag*C2(1,t)*Af*((ub-lb)*rand+lb);%eq.(6)
            end
        end
    %% Exploitation Phase (Food Exists)
    else
        if Temp>Thresold2  %hot
            for i=1:Nm
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewm(i,j)=Xfood(j)+C3(1,t)*Flag*Temp*rand*(Xfood(j)-Xm(i,j));%eq.(7)
                end
            end
            for i=1:Nf
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewf(i,j)=Xfood(j)+Flag*C3(1,t)*Temp*rand*(Xfood(j)-Xf(i,j));%eq.(7)
                end
            end
        else %cold
            if rand>0.6 %fight
                for i=1:Nm
                    for j=1:1:dim
                        FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));
                        Xnewm(i,j)=t1(1,t)*Xm(i,j) +C3(1,t)*FM*rand*(Q*Xbest_f(j)-Xm(i,j));%eq.(8)
                    end
                end
                for i=1:Nf
                    for j=1:1:dim
                        FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));
                        Xnewf(i,j)=t2(1,t)*Xf(i,j)+C3(1,t)*FF*rand*(Q*Xbest_m(j)-Xf(i,j));%eq.(9)
                    end
                end
            else%mating
                for i=1:Nm
                    for j=1:1:dim
                        Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));
                        Xnewm(i,j)=Xm(i,j) +C3(1,t)*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(10)
                    end
                end
                for i=1:Nf
                    for j=1:1:dim
                        Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));
                        Xnewf(i,j)=Xf(i,j) +C3(1,t)*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(11)
                    end
                end 
                [~, index]=sort(fitness_m);
                [~, index1]= sort(fitness_f);%排序
                Xnewm(index(end-3),:)=lb+rand*(ub-lb);
                Xnewm(index(end-2),:)=lb+rand*(ub-lb);
                Xnewm(index(end-1),:)=lb+rand*(ub-lb);
                Xnewm(index(end),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-3),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-1),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-2),:)=lb+rand*(ub-lb);
                Xnewf(index1(end),:)=lb+rand*(ub-lb);
            end
        end
    end
    %% Return back the search agents that go beyond the boundaries of the search space
    for j=1:Nm
        Flag4ub=Xnewm(j,:)>ub;
        Flag4lb=Xnewm(j,:)<lb;
        Xnewm(j,:)=(Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewm(j,:));
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
        end
    end

    %% Return back the search agents that go beyond the boundaries of the search space
    for j=1:Nf
        Flag4ub=Xnewf(j,:)>ub;
        Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewf(j,:));
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);
        end
    end

    

         for j = 1:Nm
                %% Cauchy mutation
                Temp = GaussianMutation(Xm(j,:),lb,ub);
                %% Return back the search agents that go beyond the boundaries of the search space
                Temp(Temp>ub) = ub2(Temp>ub);
                Temp(Temp<lb) = lb2(Temp<lb);
                ftemp = fobj(Temp);
                if(ftemp<fitness_m(j))
                    fitness_m(j)= ftemp;
                    Xm(j,:) = Temp;
                end
        end


        for j = 1:Nf

            Temp = GaussianMutation(Xf(j,:),lb,ub);
            %% Return back the search agents that go beyond the boundaries of the search space
            Temp(Temp>ub) = ub2(Temp>ub);
            Temp(Temp<lb) = lb2(Temp<lb);
            ftemp = fobj(Temp);
            if(ftemp<fitness_f(j))
                fitness_f(j) = ftemp;
                Xf(j,:) = Temp;
            end
        end
   
   
        [Ybest1,gbest1] = min(fitness_m);
        BestIndex1=gbest1;
            [Ybest2,gbest2] = min(fitness_f);
  BestIndex2=gbest2;
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





