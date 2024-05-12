
function [Xfood, fval,gbest_t] = so_o(N,T,lb,ub,dim,fobj)
%initial 
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
[fitnessBest_m, gbest1] = min(fitness_m);
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
Xbest_f = Xf(gbest2,:);
for t = 1:T
  Temp=exp(-((t)/T));  %eq.(4)
  Q=C1*exp(((t-T)/(T)));%eq.(5)

    % Exploration Phase (no Food)
        disp("2");
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





