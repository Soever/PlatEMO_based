%% RLSO5_1
% 雌种群so + ConvexLensImaging，雄种群RLSO


function [Xfood,fval,gbest_t] = RLSO5_3(N,T,lb,ub,dim,fobj)

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
X=lb+rand(N,dim).*(ub-lb);%eq.(1)
fitness=zeros(1,N);
for i=1:N
       fitness(i)=fobj(X(i,:));
end


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

k = N/10 ;


RF_num = 5 ;RD_num =5;strategy_num = 4 ;
q_table_m = zeros(RF_num,RD_num,strategy_num);
X_out = [];
fitness_out=[];
failure_times_m = zeros(Nm,1) ;
failure_times_f = zeros(Nf,1) ;
maxFailure_times = 10;


%% Main loop
for t = 1:T

    Temp=exp(-((t)/T));  %eq.(4)
    Q=C1(1,t)*exp(((t-T)/(T)));%eq.(5)
    
  indices = find(failure_times_m >= maxFailure_times);
  indices1 = find(failure_times_f >= maxFailure_times);
  if ~isempty(indices) || ~isempty(indices1)
      if isempty(fitness_out)
            X_out=[Xm(indices,:);Xf(indices1, :)] ;
            fitness_out=[fitness_m(indices),fitness_f(indices1)];
      else
            X_out=[X_out;Xm(indices,:);Xf(indices1, :)] ;
            fitness_out=[fitness_out,fitness_m(indices),fitness_f(indices1)];
      end
    
      [~, sortedIndex] = sort(fitness_out);
      X_out = X_out(sortedIndex, :);fitness_out = fitness_out(sortedIndex);
      
      if length(fitness_out)> N
            X_out = X_out(1:N, :);
            fitness_out = fitness_out(1:N);
      end
      
      mF_max = 0.9 ;mF_min= 0.1 ;
      mF = mF_max - (mF_max - mF_min) * (t / T);
      F = (mF - 0.1) + (0.2) * rand(); % 其中rand()生成[0,1]之间的均匀随机数
      p = round(0.05* length(fitness_out)); if p < 1 ; p=1 ;end 
        xall=  [X_out;Xm;Xf] ;
        fitnesspbest = [fitness_out,fitness_m,fitness_f];
      [~, pIndex] = sort(fitnesspbest);
      xpbest = xall(pIndex(1:p),:) ;
      
      for i = 1:length(indices)
            r =  randperm(size(xall,1), 2) ;
            x_selected = xall(r,:);
            randp = randi([1, p]) ;
            tempXm =  xpbest(randp,:)+F*(x_selected(1,:)-x_selected(2,:));
            Flag4ub=tempXm(1,:)>ub;
            Flag4lb=tempXm(1,:)<lb;
            tempXm(1,:)=(tempXm(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            y = feval(fobj,tempXm(1,:));
            fitness_m(indices(i))=y;
            Xm(indices(i),:)= tempXm;
      end
      for i = 1:length(indices1)
            r =  randperm(size(xall,1), 2) ;
            x_selected = xall(r,:);
            randp = randi([1, p]) ;
            tempXf(1,:) =  xpbest(randp,:)+F*(x_selected(1,:)-x_selected(2,:));
            Flag4ub=tempXf(1,:)>ub;
            Flag4lb=tempXf(1,:)<lb;
            tempXf(1,:)=(tempXf(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            y = feval(fobj,tempXf(1,:));
            fitness_f(indices1(i))=y;
            Xf(indices1(i),:)= tempXf;
            
    % 生成均匀分布的随机变量F
      end
  end




    state_m = get_state_knn(Xm,fitness_m,k);
    action_m = get_action(q_table_m,state_m) ;
    newXm_dec = zeros(size(Xm));

    [~, index]=sort(fitness_m);
    [~, index1]= sort(fitness_f);%排序


    for i = 1: Nm
        
        if action_m(i)==1
             newXm_dec(i,:) = ind_exploration_NoFood(Xm,fitness_m,fitness_m(i),C2(1,t),lb,ub);
        elseif action_m(i)==2
             newXm_dec(i,:) = exploit_Food(Xm(i,:),Xfood,Temp,C3(1,t));
        elseif action_m(i)==3 
            newXm_dec(i,:) = so_fight(Xm(i,:),fitness_m(i),Xbest_f,fitnessBest_f,t1(1,t),C3(1,t),Q) ;
        else
            [newXm_dec_, ~] = so_mating(Xm(i,:),Xf(i,:),fitness_m(i),fitness_f(i),C3(1,t),Q,lb,ub);
            newXm_dec(i,:) = newXm_dec_;
        end
    end
    
    if Q<Threshold
        newXf_dec = exploration_NoFood(Xf,fitness_f,C2(1,t),lb,ub);
    else
        if Temp>Thresold2
            newXf_dec = exploit_Food(Xf,Xfood,Temp,C3(1,t));
        else
            if rand > 0.6
                newXf_dec = so_fight(Xf,fitness_f,Xbest_m,fitnessBest_m,t1(1,t),C3(1,t),Q) ;
            else
                [~, newXf_dec] = so_mating(Xm,Xf,fitness_m,fitness_f,C3(1,t),Q,lb,ub);
            end
        end
    end
    for i = 0:round(Nf/10)
        newXm_dec(index(end-i),:)=lb+rand*(ub-lb);
        newXf_dec(index1(end-i),:)=lb+rand*(ub-lb);
    end
    [Xm,fitness_m,reward_m,failure_times_m] = Evaluation_reward(Xm,newXm_dec,fitness_m,lb,ub,fobj,failure_times_m);
    [Xf,fitness_f,~,failure_times_f] = Evaluation_reward(Xf,newXf_dec,fitness_f,lb,ub,fobj,failure_times_f);
    next_state_m = get_state_knn(Xm,fitness_m,k);

    for i =1:Nm
        q_table_m = updataQtable(state_m(i,:),action_m(i),reward_m(i),next_state_m(i,:),q_table_m);
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
function [newX,fitness,reward,failure_times] = Evaluation_reward(X,newX_dec,fitness,lb,ub,fobj,failure_times)
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
            failure_times(j) = 0; 
        else
            failure_times(j) = failure_times(j)+1 ; 
        
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

function [X_dec,fitness]=Evaluation_CheckBound(newX_dec,lb,ub,fobj)
    N = size(newX_dec,1);
    fitness  = zeros(1,N); 
    for i = 1:N
        Flag4ub=newX_dec(1,:)>ub;
        Flag4lb=newX_dec(1,:)<lb;
        newX_dec(1,:)=(newX_dec(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,newX_dec(1,:));
    end
    X_dec = newX_dec;

end