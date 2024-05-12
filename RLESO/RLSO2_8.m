%% RLSO2_8
% state 种群的平均适应度，种群的多样性 3*3
% action so的四个策略 （都执行择优）
% reward 
function [Xfood, fval,gbest_t,Trajectories,fitness_history, position_history] = RLSO2_8(N,T,lb,ub,dim,fobj)

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
BestIndex1=gbest1; BestIndex2=gbest2;
%% state
fitness_average_t = zeros(2,T) ;
diversity_average_t = zeros(2,T) ;
Diagonal_Length = sqrt(sum((ub2  - lb2 ).^2));

diversity_m = cal_diversity(Xm) ;
diversity_f = cal_diversity(Xf) ;

d0_m = sum(diversity_m) / (Nm*Diagonal_Length) ;
f0_m = mean(fitness_m)  ;
d0_f = sum(diversity_f) / (Nf*Diagonal_Length) ;
f0_f = mean(fitness_f)  ;
RF_num = 3 ;RD_num = 3;strategy_num = 4 ;
q_table_m = zeros(RF_num,RD_num,strategy_num);
q_table_f = zeros(RF_num,RD_num,strategy_num);



%% Main loop
for t = 1:T
    if t ==699
        a = 1 ;
    end
    disp(t);
    Temp=exp(-((t)/T));  %eq.(4)
    Q=C1(1,t)*exp(((t-T)/(T)));%eq.(5)
    Positions=[Xm;Xf];
    for i=1:size(Positions,1)
        position_history(i,t,:)=Positions(i,:);
        Trajectories(:,t)=Positions(:,1);
        fitness_history(i,t)=fobj(Positions(i,:));
    end
    k = 10*(1-2*(t/T)^2);% scaling factor eq.(19)
    
    state_m = get_state_3(Xm,fitness_m,d0_m,f0_m,Diagonal_Length);
    state_f = get_state_3(Xf,fitness_f,d0_f,f0_f,Diagonal_Length);
    action_m = get_action(q_table_m,state_m) ;
    action_f = get_action(q_table_f,state_f) ;
    if action_m==1
         newXm_dec = exploration_NoFood(Xm,fitness_m,C2(1,t),lb,ub);
    elseif action_m==2
         newXm_dec = exploit_Food(Xm,Xfood,Temp,C3(1,t));
    elseif action_m==3 
        newXm_dec = so_fight(Xm,fitness_m,Xbest_f,fitnessBest_f,t1(1,t),C3(1,t),Q) ;
    else
        [newXm_dec, ~] = so_mating(Xm,Xf,fitness_m,fitness_f,C3(1,t),Q,lb,ub);
    end
    if action_f==1
         newXf_dec = exploration_NoFood(Xf,fitness_f,C2(1,t),lb,ub);
    elseif action_f==2
         newXf_dec = exploit_Food(Xf,Xfood,Temp,C3(1,t));
    elseif action_f==3 
        newXf_dec = so_fight(Xf,fitness_f,Xbest_m,fitnessBest_m,t1(1,t),C3(1,t),Q) ;
    else
        [~, newXf_dec] = so_mating(Xm,Xf,fitness_m,fitness_f,C3(1,t),Q,lb,ub);
    end
    [Xm,fitness_m,reward_m] = Evaluation_reward(Xm,newXm_dec,fitness_m,lb,ub,fobj);
    [Xf,fitness_f,reward_f] = Evaluation_reward(Xf,newXf_dec,fitness_f,lb,ub,fobj);
    next_state_m = get_state_3(Xm,fitness_m,d0_m,f0_m,Diagonal_Length);
    next_state_f = get_state_3(Xf,fitness_f,d0_f,f0_f,Diagonal_Length);
    q_table_m = updataQtable(state_m,action_m,reward_m,next_state_m,q_table_m);
    q_table_f = updataQtable(state_f,action_f,reward_f,next_state_f,q_table_f);

    [Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f,GYbest,Xfood] = updateXbest(Xm,Xf,fitness_m,fitness_f,Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f);
    gbest_t(1,t) = GYbest ;
end
    fval = GYbest;
    showQ_table(q_table_m);
    showQ_table(q_table_f);
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
    diversity = cal_diversity(X) ;
    D = sum(diversity) / (N*DL) ;
    F = mean(fitness)  ;
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

function state = get_state_3(X,fitness,d0,f0,DL)
    [N,dim] = size(X) ;
    state = zeros(1,2) ;
    diversity = cal_diversity(X) ;
    D = sum(diversity) / (N*DL) ;
    F = mean(fitness)  ;
    RD = D/d0 ;
    RF = F/f0 ;
    if RD < 0.5
        state(1,2) = 1 ;
    elseif RD < 1
        state(1,2) = 2 ;
    else
        state(1,2) = 3;
    end
    if RF < 0.5
        state(1,1) = 1 ;
    elseif RF < 1
        state(1,1) = 2 ;
    else
        state(1,1) = 3;
    end
    
end



function Probability = softmax(x)
    % x的size为(1,n)
    % 计算每个元素的指数
    exp_x = exp(x - max(x));  % 减去max(x)增加数值稳定性
    % 计算Softmax
    Probability = exp_x / sum(exp_x);
end

function [X_dec,fitness,reward,next_state]  = act(action,action2,X_dec,X2,fitness,fitness2,N,lb,ub,fobj,d0_m,f0_m,Diagonal_Length,t,C2,C3,Q,Temp,Xfood)
    newX_dec = zeros(size(X_dec));
    if action==1
        newX_dec = exploration_NoFood(X_dec,fitness,C2,lb,ub);
    elseif action == 2
        newX_dec = exploit_Food(X_dec,Xfood,Temp,C3);
    elseif action == 3
        newX_dec = so_fight(X_dec,fitness2,Xbest,min(fitnessBest),t1,C3,Q) ;
    elseif action == 4 && action2 ==4
        [newX_dec, ~] = so_mating(X_dec,X2,fitness,fitness2,C3,Q,lb,ub);
    else
        error("action error");
    end
    fitness_new  = zeros(size(fitness)) ;
    for j=1:N
        Flag4ub=newX_dec(j,:)>ub;
        Flag4lb=newX_dec(j,:)<lb;
        newX_dec(j,:)=(newX_dec(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        fitness_new(j) = feval(fobj,newX_dec(j,:));
        
    end
    alpha = 0.5 ;
    reward  =alpha* (min(fitness) -min(fitness_new))+ (1-alpha)*(mean(fitness)-mean(fitness_new))  ;
    next_state = get_state(newX_dec,fitness_new,d0_m,f0_m,Diagonal_Length);
    X_dec=newX_dec ;
    fitness=fitness_new;
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


function diversity = cal_diversity(X_dec)
    [N,dim] = size(X_dec);
    d_pop = 0 ;
    diversity = zeros(N,1) ;
    x_mean = mean(X_dec);
    for i = 1:N
        d_ind = 0 ;
        for j = 1:dim
            d_ind = d_ind+(X_dec(i,j) - x_mean(1,j))^2; 
        end
        d_pop  = d_pop + d_ind  ;
        diversity(i,1) = d_pop ;
    end
end