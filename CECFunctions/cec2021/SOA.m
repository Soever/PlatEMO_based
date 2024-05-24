%_________________________________________________________________________%
%º£Å¸Ëã·¨             %
%_________________________________________________________________________%
function [Best_pos,Best_score,curve]=SOA(pop,Max_iter,lb,ub,dim,fobj)

fc = 2;%¿Éµ÷

if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end

%ÖÖÈº³õÊ¼»¯
X0=initialization(pop,dim,ub,lb);
X = X0;
%¼ÆËã³õÊ¼ÊÊÓ¦¶ÈÖµ
fitness = zeros(1,pop);
for i = 1:pop
   fitness(i) =  fobj(X(i,:));
end
 [fitness, index]= sort(fitness);%ÅÅÐò
GBestF = fitness(1);%È«¾Ö×îÓÅÊÊÓ¦¶ÈÖµ
%°´ÊÊÓ¦¶ÈÅÅÐò,X(1,:)´ú±í×îÓÅÎ»ÖÃ£¬X(end,:)´ú±í×î²îÎ»ÖÃ
for i = 1:pop
    X(i,:) = X0(index(i),:);
end

GBestX = X(1,:);%È«¾Ö×îÓÅÎ»ÖÃ
curve=zeros(1,Max_iter);
X_new = X;
Ms = zeros(pop,dim);
Cs = zeros(pop,dim);
Ds = zeros(pop,dim);
for t = 1: Max_iter
    
    Pbest = X(1,:);
   for i = 1:pop
        %% ¼ÆËãCs
        A = fc - (t*(fc/Max_iter));   
        Cs(i,:) =  X(i,:).*A;

        %% ¼ÆËãMs
        rd = rand(1,dim);
        B = 2*A^2*rd;
        Ms(i,:) = B.*(Pbest - X(i,:));

        %% ¼ÆËãDs
        Ds(i,:) = abs(Cs(i,:) + Ms(i,:));

        %% ¾Ö²¿ËÑË÷
        u = 1; v = 1;
        theta = rand(1,dim);
        r = u.*exp(theta*v);
        x = r.*cos(theta.*2.*pi);
        y = r.*sin(theta.*2.*pi);
        z = r.*theta;
    %% Î»ÖÃ¸üÐÂ
        X_new(i,:) = x.*y.*z.*Ds(i,:) + Pbest;
    end
   %±ß½ç¿ØÖÆ
   for j = 1:pop
       for a = 1: dim
           if(X_new(j,a)>ub)
               X_new(j,a) =ub(a);
           end
           if(X_new(j,a)<lb)
               X_new(j,a) =lb(a);
           end
       end
   end 
   %¸üÐÂÎ»ÖÃ
   for j=1:pop
    fitness_new(j) = fobj(X_new(j,:));
   end
   for j = 1:pop
    if(fitness_new(j) < GBestF)
        GBestF = fitness_new(j);
        GBestX = X_new(j,:);   
    end
   end
    X = X_new;
    fitness = fitness_new;
    %ÅÅÐò¸üÐÂ
   [fitness, index]= sort(fitness);%ÅÅÐò
   for j = 1:pop
      X(j,:) = X(index(j),:);
   end
   curve(t) = GBestF;
end
Best_pos = GBestX;
Best_score = curve(end);
end



