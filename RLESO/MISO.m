
function [Xfood, fval,gbest_t ] = MISO(N,T,lb,ub,dim,fobj)
%initial 
history = zeros(1,T);
vec_flag=[1,-1];
Threshold=0.25;
Thresold2= 0.6;
C1=0.5;
C2=.05;
C3=2;
X=lb+rand(N,dim)*(ub-lb);
for i=1:N
    fitness(i)=feval(fobj,X(i,:));   
end
[GYbest, gbest] = min(fitness);
Xfood = X(gbest,:);
%Diving the swarm into two equal groups males and females
Nm=round(N/2);%eq.(2&3)
Nf=N-Nm;
Xm=X(1:Nm,:);
Xf=X(Nm+1:N,:);
fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:N);
X_out = [];
fitness_out=[];
failure_times_m = zeros(Nm,1) ;
failure_times_f = zeros(Nf,1) ;
maxFailure_times = 5 ;
[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
for t = 1:T
  Temp=exp(-((t)/T));  %eq.(4)
  Q=C1*exp(((t-T)/(T)));%eq.(5)

   

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
            disp("maxFailure");
      end
      F = (mF - 0.1) + (0.2) * rand();
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
      end
  end
    if Q>1        Q=1;    end
    % Exploration Phase (no Food)
if Q<Threshold
    for i=1:Nm
        for j=1:1:dim
            rand_leader_index = floor(Nm*rand()+1);
            X_randm = Xm(rand_leader_index, :);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));%eq.(7)
            Xnewm(i,j)=X_randm(j)+Flag*C2*Am*((ub-lb)*rand+lb);%eq.(6)
        end
    end
    for i=1:Nf
        for j=1:1:dim
            rand_leader_index = floor(Nf*rand()+1);
            X_randf = Xf(rand_leader_index, :);
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));%eq.(9)
            Xnewf(i,j)=X_randf(j)+Flag*C2*Af*((ub-lb)*rand+lb);%eq.(8)
        end
    end
else %Exploitation Phase (Food Exists)
    if Temp>Thresold2  %hot
        for i=1:Nm
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewm(i,j)=Xfood(j)+C3*Flag*Temp*rand*(Xfood(j)-Xm(i,j));%eq.(10)
            end
        end
        for i=1:Nf
            flag_index = floor(2*rand()+1);
            Flag=vec_flag(flag_index);
            for j=1:1:dim
                Xnewf(i,j)=Xfood(j)+Flag*C3*Temp*rand*(Xfood(j)-Xf(i,j));%eq.(10)
            end
        end
    else %cold
        if rand>0.6 %fight
            for i=1:Nm
                for j=1:1:dim
                    FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));%eq.(13)
                    Xnewm(i,j)=Xm(i,j) +C3*FM*rand*(Q*Xbest_f(j)-Xm(i,j));%eq.(11)
                end
            end
            for i=1:Nf
                for j=1:1:dim
                    FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));%eq.(14)
                    Xnewf(i,j)=Xf(i,j)+C3*FF*rand*(Q*Xbest_m(j)-Xf(i,j));%eq.(12)
                end
            end
        else%mating
            for i=1:Nm
                for j=1:1:dim
                    Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));%eq.(17)
                    Xnewm(i,j)=Xm(i,j) +C3*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(15
                end
            end
            for i=1:Nf
                for j=1:1:dim
                    Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));%eq.(18)
                    Xnewf(i,j)=Xf(i,j) +C3*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(16)
                end
            end
            flag_index = floor(2*rand()+1);
            egg=vec_flag(flag_index);
            if egg==1;
                [GYworst, gworst] = max(fitness_m);
                Xnewm(gworst,:)=lb+rand*(ub-lb);%eq.(19)
                [GYworst, gworst] = max(fitness_f);
                Xnewf(gworst,:)=lb+rand*(ub-lb);%eq.(20)
            end
        end
    end
end

 hm = zeros(1,Nm);hf = zeros(1,Nm);
    for j=1:Nm
         Flag4ub=Xnewm(j,:)>ub;
         Flag4lb=Xnewm(j,:)<lb;
        Xnewm(j,:)=(Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewm(j,:));
        hm(j) = y ; 
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
            failure_times_m(j) = 0; 
        else
            failure_times_m(j) = failure_times_m(j)+1 ; 
        end
    end
    
    [Ybest1,gbest1] = min(fitness_m);
    
    for j=1:Nf
         Flag4ub=Xnewf(j,:)>ub;
         Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewf(j,:));
         hf(j) = y ;
        if y<fitness_f(j)
             
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);
            failure_times_f(j) = 0; 
        else
            failure_times_f(j) = failure_times_f(j)+1 ; 
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





