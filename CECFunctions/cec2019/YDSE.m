%_______________________________________________________________________________________%
%Young's Double-Slit Experoment (YDSE) Optimizer
%paper: Abdel-Basset et al., "Young's double-slit experiment Optimizer for global and constraint optimization problems".
%published in Computer methods in applied mechanics and engineering,2022
%F_obj--------objective function
%NP-------population size
%Dim-------dimension of problems
%LB-----the lower limit of the vaiables
%UB-----the upper limit of the vaiables
%Max_iter------the maximum number of function evaluations
%cgcurve---the record of the convergence curves
%Best_Fringe_fitness--the optimal Fitness value
%BestFringe------the optimal solution
%Developed in MATLAB R2019a by Doaa El-Shahat & Mohamed Abdel-Basset                                                   %
%_______________________________________________________________________________________%
function [Best_Fringe_fitness,BestFringe,cgcurve]=YDSE(NP,Max_iter,LB,UB,Dim,F_obj)
%--------------------------------------------------------------------------
%         YDSE parameters
%--------------------------------------------------------------------------
L=1; % L is the distance between the barrier and the screen is often in the order of 1 m
d=5*10^-3; %distance between two slits FS and SS is often a fraction of millimeter
I=0.01; % distance between light source and barrier
Lambda=5*10^-6; % is the wavelength is a fraction of micrometer
Delta=0.38; % is a constant
%--------------------------------------------------------------------------
% Phase 1: initialize a source of Monochromatic light waves
%--------------------------------------------------------------------------
S=initialization(NP,Dim,UB,LB); % Eq. (6) Initialize population  monochromatic light is projected S using Eq.(4) and Eq.(5)
%--------------------------------------------------------------------------
%Phase 2: Huygens' principle to generate the waves from the double-slits FS and SS
%--------------------------------------------------------------------------
FS=zeros(NP,Dim);
SS=zeros(NP,Dim);
S_Mean=mean(S);
for i=1:NP
    FS(i,:)= S(i,:)+I.*(2*rand-1).*(S_Mean-S(i,:)) ; %Eq. (9)
    SS(i,:)= S(i,:)-I.*(2*rand-1).*(S_Mean-S(i,:)) ; %Eq. (10)
end
%--------------------------------------------------------------------------
% Phase 3: travelling Two waves outgoing from the double slits : path difference
%--------------------------------------------------------------------------
BestFringe=zeros(1,Dim);   % the best-obtained so far solution
Best_Fringe_fitness = inf;      % the optimal Fitness value
Fringes_fitness=zeros(1,NP);
X=zeros(NP,Dim);
for i=1:NP
    if i==1   %m=i-1   then m=0   central fringe
        Delta_L=0;% AS THE TWO WAVES HAS THE SAME LEGTH  .. no path difference
        for j=1:Dim
            X(i,j)=(FS(i,j)+SS(i,j))/2+ Delta_L; %Eq.(12)
        end
    elseif mod(i,2)==0 %odd values of m=1, 3, 5, 7,...
        m=i-1;
        Delta_L=(2*m+1)*Lambda/2; % compute path difference for the odd m number
        for j=1:Dim
            X(i,j)= (FS(i,j)+SS(i,j))/2+Delta_L;%Eq.(12)
        end
    else            %even values of m=2,4,6,... compute path difference for the even m number
        m=i-1; %order number
        Delta_L=(m)*Lambda;
        for j=1:Dim
            X(i,j)= (FS(i,j)+SS(i,j))/2+Delta_L;%Eq.(12)
        end
    end
end
for i=1:NP
    F_UB=X(i,:)>UB;
    F_LB=X(i,:)<LB;
    X(i,:)=(X(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
end
for i=1:NP
    Fringes_fitness(1,i) = F_obj(X(i,:));
end
%rank the solutions according to the fitness function from best to worst
[Fringes_fitness,sorted_indices]=sort(Fringes_fitness);
temp_X=X;
X=temp_X(sorted_indices,:);
Best_Fringe_fitness=Fringes_fitness(1,1);
BestFringe = X(1,:);
%--------------------------------------------------------------------------
%Phase 4: Interference pattern on the projection screen . bright and dark fringes
%--------------------------------------------------------------------------
Int_max0=10^-20;
iter=1;
X_New=zeros(NP,Dim);
while iter<Max_iter+1
    cgcurve(iter)=Best_Fringe_fitness;
    %rank the solutions according to the fitness function sort fitness by indices
    [Fringes_fitness,sorted_indices]=sort(Fringes_fitness);
    temp_X=X;
    X=temp_X(sorted_indices,:);
    q=iter/Max_iter;   %Eq.(26)
    Int_max=Int_max0.*q;   %Eq.(25)
    for  i=1:NP
        a=iter.^(2*rand()-1); %Eq.(20)
        H=2.*rand(1,Dim)-1;
        z=a./H;   %Eq.(19)
        r1=rand;
        %------------exploitation phase----------------
        if i==1 %zero m value
            beta=q.*cosh(pi./iter); %Eq.(15)
            A_bright=2./(1+sqrt(abs(1-beta.^2))); %Eq.(14)
            A=1:2:NP;
            randomIndex = randi(length(A));%select random bright fringe
            X_New(i,:)=BestFringe+(Int_max.*A_bright.*X(i,:)-r1.*z.*X(A(randomIndex),:)); % %Eq. (24)
        elseif  mod(i,2)==1 %even m values
            beta=q.*cosh(pi./iter);%Eq.(15)
            A_bright=2./(1+sqrt(abs(1-beta.^2)));  %Eq.(14)
            m=i-1;
            y_bright=Lambda*L*m/d; %Eq.(4)
            Int_bright=Int_max.*cos((pi*d)/(Lambda*L)*y_bright).^2;   %Eq.(23)  %update light intensity using Eq.(18)
            s=randperm(NP,2);
            g=2*rand-1;
            Y=(X(s(2),:)-X(s(1),:)); %select two random fringes Eq.(22)
            X_New(i,:)=X(i,:)-((1-g).*A_bright.*Int_bright.*X(i,:)+g.*Y);   %Eq.(21)
            
            %------------exploration-------------------
        else % odd m values
            A_dark=Delta*atanh(-(q)+1); %Eq.(16)
            m=i-1;
            y_dark=Lambda*L*(m+0.5)/d; %Eq.(5)
            Int_dark=Int_max.*cos((pi*d)/(Lambda*L)*y_dark).^2;%Eq.(18)
            X_New(i,:)=X(i,:)-(r1.*A_dark.*Int_dark.*X(i,:)-z.*BestFringe); %Eq. (17)
        end
        %check bounds of new solutions
        F_UB=X_New(i,:)>UB;
        F_LB=X_New(i,:)<LB;
        X_New(i,:)=(X_New(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
    end
    %Generate the next population
    X_New_Fitness=zeros(1,size(X,1));
    for i=1:NP
        X_New_Fitness(1,i) = F_obj(X_New(i,:));
        if X_New_Fitness(1,i)<Fringes_fitness(1,i)
            X(i,:)=X_New(i,:);
            Fringes_fitness(1,i)=X_New_Fitness(1,i);
            if Fringes_fitness(1,i) < Best_Fringe_fitness
                Best_Fringe_fitness=Fringes_fitness(1,i);
                BestFringe = X(i,:);
            end
        end
    end
    %     if mod(iter,100)==0
    %         disp(['Iteration ' num2str(iter) ': Best Cost = ' num2str(Best_Fringe_fitness,20)]);
    %     end
    iter=iter+1;
end
end