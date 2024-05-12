
function [Xfood, fval,gbest_t] = so_3(N,T,lb,ub,dim,fobj)
%initial 
vec_flag=[1,-1];
Threshold=0.25;
Thresold2= 0.6;
C1=0.5;
C2=.05;
C3=2;
% C1=0.5*ones(1,T);
% C2=0.05*ones(1,T);
% C3=2*ones(1,T);


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
[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
for t = 1:T
  Temp=exp(-((t)/T));  %eq.(4)
  Q=C1*exp(((t-T)/(T)));%eq.(5)
    if Q>1        Q=1;    end
    disp("3");
    %[Xnewm, Xnewf] =so_mating(Xm,Xf,fitness_m,fitness_f,C3,Q,lb,ub);
     Xnewm = exploit_Food(Xm,Xfood,Temp,C3);
     Xnewf = exploit_Food(Xm,Xfood,Temp,C3);
    %[Xnewm, Xnewf] = so_mating(Xm,Xf,fitness_m,fitness_f,C3,Q,lb,ub);
    % Xnewm = so_fight(Xm,fitness_m,Xbest_f,fitnessBest_f,1,C3,Q) ;
     % Xnewf = so_fight(Xf,fitness_f,Xbest_m,fitnessBest_m,1,C3,Q) ;
    % for i=1:Nm
    %         flag_index = floor(2*rand()+1);
    %         Flag=vec_flag(flag_index);
    %         for j=1:1:dim
    %             Xnewm(i,j)=Xfood(j)+C3*Flag*Temp*rand*(Xfood(j)-Xm(i,j));%eq.(10)
    %         end
    %     end
    %     for i=1:Nf
    %         flag_index = floor(2*rand()+1);
    %         Flag=vec_flag(flag_index);
    %         for j=1:1:dim
    %             Xnewf(i,j)=Xfood(j)+Flag*C3*Temp*rand*(Xfood(j)-Xf(i,j));%eq.(10)
    %         end
    %     end

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
    
    [Ybest1,gbest1] = min(fitness_m);
    
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





