%% RLSO2_8
% state >mean(fitness),kean(decs)
% action so的四个策略 （都执行择优）
% reward 
function [Xfood,fval,gbest_t] = RLSO2_19(N,T,lb,ub,dim,fobj)

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
       fitness(i)=fobj(X(i,:));
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
% BestIndex1=gbest1; BestIndex2=gbest2;
%% state
% fitness_average_t = zeros(2,T) ;
% diversity_average_t = zeros(2,T) ;
% Diagonal_Length = sqrt(sum((ub2  - lb2 ).^2));
k = N/10 ;
% d0_m = get_neighbor_diversity(Xm,k,Diagonal_Length) ;
% d0_f = get_neighbor_diversity(Xf,k,Diagonal_Length) ;
% f0_m = fitness_m ;
% f0_f = fitness_f ;

RF_num = 5 ;RD_num =5;strategy_num = 4 ;
q_table_m = zeros(RF_num,RD_num,strategy_num);
q_table_f = zeros(RF_num,RD_num,strategy_num);



%% Main loop
for t = 1:T

    Temp=exp(-((t)/T));  %eq.(4)
    Q=C1(1,t)*exp(((t-T)/(T)));%eq.(5)
    Positions=[Xm;Xf];
    for i=1:size(Positions,1)
        position_history(i,t,:)=Positions(i,:);
        Trajectories(:,t)=Positions(:,1);
        fitness_history(i,t)=fobj(Positions(i,:));
        
    end

    state_m = get_state_knn(Xm,fitness_m,k);
    state_f = get_state_knn(Xf,fitness_f,k);
    action_m = get_action(q_table_m,state_m) ;
    action_f = get_action(q_table_f,state_f) ;
    
    newXm_dec = zeros(size(Xm));
    newXf_dec = zeros(size(Xf));

     [~, index]=sort(fitness_m);
     [~, index1]= sort(fitness_f);%排序
    
    kk = 10*(1-2*(t/T)^2);
    TempXm = ConvexLensImaging(kk,Xbest_m,ub,ub1,lb,lb1);
    fitTemp = fobj(TempXm);
    if(fitTemp<GYbest)
        fitnessBest_m=fitTemp ;
        Xbest_m = TempXm;
        Xm(index(1),:) = TempXm;
    end
    
    TempXf =ConvexLensImaging(kk,Xbest_f,ub,ub1,lb,lb1);

    fitTemp = fobj(TempXf);
    if(fitTemp<GYbest)
        fitnessBest_f=fitTemp ;
        Xbest_f = TempXf;
         Xf(index1(1),:) = TempXf;
    end
    for i = 1: Nm
        if action_m(i)==1
             newXm_dec(i,:) = exploration_NoFood(Xm(i,:),fitness_m(i),C2(1,t),lb,ub);
        elseif action_m(i)==2
            Temp = 4*rand()*(1-rand());
             newXm_dec(i,:) = exploit_Food(Xm(i,:),Xfood,Temp,C3(1,t));
        elseif action_m(i)==3 
            Q = 4*rand()*(1-rand()) ;
            newXm_dec(i,:) = so_fight(Xm(i,:),fitness_m(i),Xbest_f,fitnessBest_f,t1(1,t),C3(1,t),Q) ;
        else
            Q = 4*rand()*(1-rand()) ;
            [newXm_dec_, ~] = so_mating(Xm(i,:),Xf(i,:),fitness_m(i),fitness_f(i),C3(1,t),Q,lb,ub);
            newXm_dec(i,:) = newXm_dec_;
        end
    end
    for i = 1: Nf
        if action_f(i)==1
             newXf_dec(i,:) = exploration_NoFood(Xf(i,:),fitness_f(i),C2(1,t),lb,ub);
        elseif action_f(i)==2
             newXf_dec(i,:) = exploit_Food(Xf(i,:),Xfood,Temp,C3(1,t));
        elseif action_f(i)==3 
            newXf_dec(i,:) = so_fight(Xf(i,:),fitness_f(i),Xbest_m,fitnessBest_m,t1(1,t),C3(1,t),Q) ;
        else
            [~, newXf_dec(i,:)] = so_mating(Xm(i,:),Xf(i,:),fitness_m(i),fitness_f(i),C3(1,t),Q,lb,ub);
        end
    end
     [~, index]=sort(fitness_m);
     [~, index1]= sort(fitness_f);%排序
     
     for i = 0:round(Nm/10)
        newXm_dec(index(end-i),:)=lb+rand*(ub-lb);
        newXf_dec(index1(end-i),:)=lb+rand*(ub-lb);
     end
     
    [Xm,fitness_m,reward_m] = Evaluation_reward(Xm,newXm_dec,fitness_m,lb,ub,fobj);
    [Xf,fitness_f,reward_f] = Evaluation_reward(Xf,newXf_dec,fitness_f,lb,ub,fobj);
    next_state_m = get_state_knn(Xm,fitness_m,k);
    next_state_f = get_state_knn(Xf,fitness_f,k);
    for i =1:Nm-round(Nm/10)
        jm = index(i);jf = index1(i);
        q_table_m = updataQtable(state_m(jm,:),action_m(jm),reward_m(jm),next_state_m(jm,:),q_table_m);
        q_table_f = updataQtable(state_f(jf,:),action_f(jf),reward_f(jf),next_state_f(jf,:),q_table_f);
    end
    
    

    [Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f,GYbest,Xfood] = updateXbest(Xm,Xf,fitness_m,fitness_f,Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f);
    gbest_t(1,t) = GYbest ;
end
    fval = GYbest;
    % showQ_table(q_table_m);
    % showQ_table(q_table_f);
    
end

function q = updataQtable(s,a,r,s_next,q)
    actions = q(s_next(1),s_next(2),:);
    [q_target_value,index] = max(actions);
    q(s(1),s(2),a) =q(s(1),s(2),a)+0.1*(r+0.9*q_target_value-q(s(1),s(2),a)) ;
end

function action = get_action(q_table,state)
    N = size(state,1);
    action = zeros(N,1);
    for i  = 1:N
        actions_value = q_table(state(i,1),state(i,2),:);
        Probability = softmax(actions_value);
        action(i,1)= randsample(1:length(Probability), 1, true, Probability );
    end
end



function state = get_state_knn(X,fitness,k)
    [N,~] = size(X) ;
    state = zeros(N,2) ;
    D =get_neighbor_diversity(X,k);
    ub2 = max(X);
    lb2 = min(X);
    DL = sqrt(sum((ub2  - lb2 ).^2));
    F = fitness ;
    popD  = sum(cal_diversity(X)) / (N*(DL+1e-160)) ;
    popF  = mean(F) ;
    % 种群丰富度可能=0 

    
    if popD  == 0 
       RDs = zeros(size(D));
    else
       RDs = D./popD;
    end
    
    % 种群适应度可能出现全为0
    if popF == 0 
        RFs = zeros(size(F'));
    else
        RFs = F'./popF ;
    end
    RDs = mapminmax(RDs',0,1);
    RFs = mapminmax(RFs',0,1);
    state(RDs < 0.2, 2) = 1;
    state(RDs >= 0.2 & RDs < 0.4, 2) = 2;
    state(RDs >= 0.4 & RDs < 0.6, 2) = 3;
    state(RDs >= 0.6 & RDs < 0.8, 2) = 4;
    state(RDs >= 0.8, 2) = 5;
    
   
    state(RFs < 0.2, 1) = 1;
    state(RFs >= 0.2 & RFs < 0.4, 1) = 2;
    state(RFs >= 0.4 & RFs < 0.6, 1) = 3;
    state(RFs >= 0.6 & RFs <= 0.8, 1) = 4;
    state(RFs > 0.8, 1) = 5;
end

function neighbor_diversities =get_neighbor_diversity(X,k)
    [N,dim] = size(X) ;
%     nearest_indices_all = zeros(N, k+1);
    [nearest_indices_all,~] = knnsearch(X,X,'k',k+1);
    % 遍历每个个体
%     for i = 1:N
%         % 使用 knnsearch 函数找到每个个体最近的 k 个个体的索引
%         nearest_indices = knnsearch(X, X(i, :), 'K', k+1); % +1 是因为包含了自身#问题1：knnsearch可以直接获取每个个体对应的近邻，不需要放到循环里
%         % 将最近的 k 个个体的索引存储到矩阵中
%         nearest_indices_all(i, :) = nearest_indices;
%     end
    neighbor_diversities = zeros(N,1);
    for i = 1:N
        neighbor_dec = X(nearest_indices_all(i,:),:) ;
        ub2 = max(neighbor_dec);
        lb2 = min(neighbor_dec);
        Diagonal_Length = sqrt(sum((ub2  - lb2 ).^2));
        neighbor_diversities(i) =sum(cal_diversity(neighbor_dec))/((k+1)*(Diagonal_Length+1e-160)) ;%要取独立的DL
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
function [newX,fitness,reward] = Evaluation_reward(X,newX_dec,fitness,lb,ub,fobj)
    N = size(X,1);
    fitness_new = zeros(size(fitness));
    fitness_old = fitness ;
    reward = -1*ones(N,1) ;
    for j=1:N
        Flag4ub=newX_dec(j,:)>ub;
        Flag4lb=newX_dec(j,:)<lb;
        newX_dec(j,:)=(newX_dec(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
       fitness_new(j)=fobj(newX_dec(j,:));
   
%         fitness_new(j) = feval(fobj,newX_dec(j,:));
        % 择优
        if fitness_new(j)  < fitness(j)
            X(j,:) = newX_dec(j,:) ;
            fitness(j) = fitness_new(j);
            reward(j) = 1 ;
        end
    end
    % reward = fitness_old-fitness_new ;
    
    
    
    newX = X;
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
        diversity(i,1) = sqrt(d_pop) ;
    end
end