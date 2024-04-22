
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
%% state
fitness_average_t = zeros(2,T) ;
diversity_average_t = zeros(2,T) ;
Diagonal_Length = sqrt(sum((ub2  - lb2 ).^2));
d0_m = sum(std(Xm,1)) / (Nm*Diagonal_Length) ;
f0_m = mean(fitness_m) / Nm ;
d0_f = sum(std(Xf,1)) / (Nf*Diagonal_Length) ;
f0_f = mean(fitness_f) / Nf ;
RF_num = 5 ;RD_num = 5 ;strategy_num = 3 ;
q_table_m = zeros(RF_num,RD_num,strategy_num);
q_table_f = zeros(RF_num,RD_num,strategy_num);
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
    state_m = get_state(Xm,fitness_m,d0_m,f0_m,Diagonal_Length);
    state_f = get_state(Xf,fitness_f,d0_f,f0_f,Diagonal_Length);
    action_m = get_action(q_table_m,state_m) ;
    action_f = get_action(q_table_f,state_f) ;
    [Xm,fitness_m,reward_m,next_state_m]= act(action_m,Xm,fitness_m,Nm,lb,ub,fobj,d0_m,f0_m,Diagonal_Length,t);
    [Xf,fitness_f,reward_f,next_state_f]= act(action_f,Xf,fitness_f,Nf,lb,ub,fobj,d0_f,f0_f,Diagonal_Length,t);
    q_table_m = updataQtable(state_m,action_m,reward_m,next_state_m,q_table_m);
    q_table_f = updataQtable(state_f,action_f,reward_f,next_state_f,q_table_f);

    [Xbest_m,Xbest_f,fitness_m,fitness_f,GYbest,Xfood] = updateXbest(Xm,Xf,fitness_m,fitness_f,Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f);
    gbest_t(1,t) = GYbest ;
end
    fval = GYbest;
end

function q = updataQtable(s,a,r,s_next,q)
    actions = q(s_next(1),s_next(2),:);
    [q_target_value,index] = max(actions);
    q(s(1),s(2),a) =q(s(1),s(2),a)+0.1*(r+0.9*q_target_value-q(s(1),s(2),a)) ;
end

function action = get_action(q_table,state)
    actions = q_table(state(1),state(2),:);
    Probability = softmax(actions);
    action= randsample(1:length(Probability), 1, true, Probability );
end
function state = get_state(X,fitness,d0,f0,DL)
    [N,dim] = size(X) ;
    state = zeros(1,2) ;
    D = sum(std(X,1)) / (N*DL) ;
    F = mean(fitness) / N ;
    RD = D/d0 ;
    RF = F/f0 ;
    
    if RD < 0.3
        state(1,2) = 1 ;
    elseif RD < 0.6
        state(1,2) = 2 ;
    elseif RD < 1
        state(1,2) = 3 ;
    elseif RD < 1.3
        state(1,2) = 4 ;
    else
        state(1,2) = 5 ;
    end
    if RF < 0.3
        state(1,1) = 1 ;
    elseif RF < 0.5
        state(1,1) = 2 ;
    elseif RF < 0.7
        state(1,1) = 3 ;
    elseif RF <= 1
        state(1,1) = 4 ;
    else
        state(1,1) = 5 ;
    end
    

    
end
function Probability = softmax(x)
    % x的size为(1,n)
    % 计算每个元素的指数
    exp_x = exp(x - max(x));  % 减去max(x)增加数值稳定性
    % 计算Softmax
    Probability = exp_x / sum(exp_x);
end

function [X_dec,fitness,reward,next_state]  = act(action,X_dec,fitness,N,lb,ub,fobj,d0,f0,Diagonal_Length,t)
    newX_dec = zeros(size(X_dec));
    if action==1
        newX_dec = CauchyMutation(X_dec);
    elseif action == 2
        newX_dec = GaussianMutation(X_dec,lb,ub);
    elseif action == 3
        newX_dec = tMutation(X_dec,t);
    else
        error("action error");
    end
    fitness_new  = zeros(size(fitness)) ;
    for j=1:N
        Flag4ub=newX_dec(j,:)>ub;
        Flag4lb=newX_dec(j,:)<lb;
        newX_dec(j,:)=(newX_dec(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness_new(j) = feval(fobj,newX_dec(j,:));
        if fitness_new(j) <fitness(j)
            fitness(j) = fitness_new(j);
            X_dec(j,:) = newX_dec(j,:) ;
            reward(j) = 1 ;
        else
            reward(j) = -1 ;
        end
    end
    reward = sum(reward) ;
    next_state = get_state(newX_dec,fitness_new,d0,f0,Diagonal_Length);
     
end
function [X,fitness,reward] = Evaluation_reward(X,newX_dec,fitness,lb,ub,fobj)
    N = size(X,1);
    reward = zeros(1,N);
    for j=1:N
        Flag4ub=newX_dec(j,:)>ub;
        Flag4lb=newX_dec(j,:)<lb;
        newX_dec(j,:)=(newX_dec(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,newX_dec(j,:));
        if y<fitness(j)
            fitness(j)=y;
            X(j,:)= newX_dec(j,:);
            reward(1,j) = 1 ;
        else
            reward(1,j) = -1 ;
        end

    end
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