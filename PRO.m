clc;
%%
   Func_id = 1;          % CEC2017 1~30 
   nPop = 50;            % Population Size
   N = 10;               % number of Decision variables 
   MaxFEs = N * 10000;   % Maximum number of function evaluations
   NumofExper = 1;       % Number of test
  %Benchmark = 3;        % 1:Classic 2:CEC2005 3:CEC2014 
   LB=-100;%lb;          % Lower Bound
   UB=100;%ub;           % Upper Bound
% 
global initial_flag
initial_flag = 0;
%% 
Function_name=['F' num2str(Func_id)];
%========== CEC2017 ==========
%fhd=str2func('cec14_func');
fhd=str2func('cec17_func');
CostFunction=Func_id;
%============================= 
LB = LB.*ones(1,N);       
UB = UB.*ones(1,N);       

% Empty Solution Structure
empty_Solution.Position=[];
empty_Solution.Cost=[];
 
Population=repmat(empty_Solution,nPop,1);
SumBestCostPRO_=zeros(MaxFEs,1);
BestSolCostPRO= []; %zeros(MaxFEs,1);

%===================================================

for ii=1:NumofExper
    
  rand('state',sum(100*clock));
  initial_flag = 0; % should set the flag to 0 for each run, each function
  
   % Create Initial Population
    for i=1:nPop
        Population(i).Position= LB+rand(1,N).*(UB-LB); %LB+rand(1,N).*(UB-LB);    
        %Population(i).Cost = feval(fhd,Population(i).Position',CostFunction) -  (CostFunction*100); % CEC2014 F(X) - F(X*)
        %Population(i).Cost = feval(fhd,Population(i).Position',CostFunction) ; % CEC2017 
        Population(i).Cost = Sphere(Population(i).Position); % Test Func
    end  
    [BestCostPRO_,BestSolCostPRO(ii)]=PRO_v1(N,MaxFEs,LB,UB,Population,nPop,CostFunction); 
end
function [BestCosts,BestSolCost]=PRO_v1(N,MaxFEs,LB,UB,population,pop_size,ObjFunc_ID)

%% Initialization 
rand('state',sum(100*clock));
% --- Problem Definition ---
                                  
%Obj_Func = @ YourObjFunc;         %  Objective Function
fhd = str2func('cec17_func');
nVar = N;                          %  Number of Decision Variables
%LB =  LB .* ones(1,nVar);          %  Variables Lower Bound
%UB =  UB .* ones(1,nVar);          %  Variables Upper Bound

% --- PRO Parameters ---
RR = 0.7;                          %  Reinforcement Rate (RR)
%MaxFEs = MaxFEs;                  %  Maximum Number of Function Evaluations
nPop = pop_size;                   %  Population Size
FEs = 0;                           %  Function Evaluations counter

% --- Empty Structure for Individuals ---
empty_individual.Behaviors=[];
empty_individual.response=[];
empty_individual.Schedule=[];
%empty_individual.Fy=[];
%empty_individual.Gx=[];

% --- Initialize Population Array ---
pop = repmat(empty_individual, nPop, 1);
% --- Initialize Best Solution ---
BestSol.response = inf;
% --- Initialize Population --- 
for i=1:nPop
    pop(i).Behaviors = population(i).Position;
    pop(i).Schedule = unifrnd(0.9,1,1,N);
    pop(i).response =  population(i).Cost; %feval('cec14_func',pop(i).Behaviors',CostFunction)  - (CostFunction*100);
end
% --- Sort pop ---
[~,SorteIndx] = sort([pop.response]);
pop = pop(SorteIndx);
% --- Set the Best Solution ---
BestSol = pop(1);

% --- Initialize Best Cost Record ---
BestCosts = zeros(MaxFEs,1);
BestCosts(1) = BestSol.response;
[~, sortedIndx] = sort([pop.response]);

ResetZero = zeros(1,N);

%% --- PRO Main Loop ---
 while FEs < MaxFEs 
    for i=1:nPop  % For all Learners      
       tempBehav = pop(i);% empty_individual;      
       k = nPop;
       if i < nPop
         k = sortedIndx(randi([i+1 nPop]));
       end    
       %% --- Determine Behaviors of the ith learner based on Scheduler. -----------  
    
       % According to Eq.(1) & Eq.(2)
       
       Tau = (FEs/MaxFEs);                          %  Time parameter
       %Selection_rate = Tau^0.5; 
       Selection_rate = exp(-(1-Tau)); %******** 
       %Selection_rate = exp(-(Tau))^2;
       [~,Candid_Behavs] = sort(pop(i).Schedule(1:N),'descend');
       % --Select Landa number of Behaviors with highest priority in Schedule i.-- 
       Landa = ceil(N*rand*Selection_rate);
       Selected_behaviors = Candid_Behavs(1:Landa);%
                
       %% --- Stimulate the selected Behaviors of the ith learner to get response.---  
       % According to Eq.(3), Eq.(4), and Eq.(5)
     
       if rand <  0.5 %(0.1 + 0.9 * (1-Tau))
          Stimulation = ResetZero;
          Stimulation(Selected_behaviors) = ( BestSol.Behaviors(Selected_behaviors) - pop(i).Behaviors(Selected_behaviors));    
       else
          Stimulation = ResetZero;
          Stimulation(Selected_behaviors) = ( pop(i).Behaviors(Selected_behaviors) - pop(k).Behaviors(Selected_behaviors));            
       end
       
       % ---- Calculate Stimulation Factor (SF) ------
       %SF = rand * ( exp(-1 * mean( abs(BestSol.Behaviors - pop(i).Behaviors)/max(abs(pop(1).Behaviors - pop(nPop).Behaviors)))));
       %SF = rand * ( exp(-1 * mean( abs(BestSol.Behaviors(Selected_behaviors) - pop(i).Schedule(Selected_behaviors))/max(abs(pop(i).Schedule(Selected_behaviors) - pop(nPop).Schedule(Selected_behaviors))))));
       SF = Tau + rand * (mean((pop(i).Schedule(Selected_behaviors) )/max(abs(pop(i).Schedule)))); %(exp(-(1-FEs/MaxFEs))^2 ) ;
       
       tempBehav.Behaviors(Selected_behaviors) = pop(i).Behaviors(Selected_behaviors) + SF .* Stimulation(Selected_behaviors);    
       
       % ------------  Bound constraints control ------------------- 
       % 
       [~,underLB] = find(tempBehav.Behaviors < LB);
       [~,uperUB] = find(tempBehav.Behaviors > UB);
       if ~isempty(underLB)
         tempBehav.Behaviors(underLB) =  LB(underLB) + rand(1,size(underLB,2)).*((UB(underLB) -  LB(underLB))./1); 
       end
       if ~isempty(uperUB)
         tempBehav.Behaviors(uperUB) =  LB(uperUB) + rand(1,size(uperUB,2)).*((UB(uperUB) -  LB(uperUB))./1); 
       end
      
       % ------ Evaluate the ith learner Response -------------------
       %tempBehav.response = feval('cec14_func',tempBehav.Behaviors',ObjFunc_ID)  - (ObjFunc_ID*100);
       
       %tempBehav.response = feval(fhd,tempBehav.Behaviors',ObjFunc_ID); %CEC2017
       tempBehav.response = Sphere(tempBehav.Behaviors'); % Test Func
     
       FEs = FEs + 1;
       
       % ----- Apply Positive or Negative Reinforcement according to the response.
      
       % According to  Eq.(6)& Eq.(7)
      
       if tempBehav.response<pop(i).response
            % Positive Reinforcement 
            tempBehav.Schedule(Selected_behaviors) = pop(i).Schedule(Selected_behaviors) + pop(i).Schedule(Selected_behaviors) * (RR/2);           
            
            % accept new Solution
            pop(i) = tempBehav;
            
            % Update the best Solution
            if pop(i).response < BestSol.response
                BestSol = pop(i);  
            end
       else
           % Negative Reinforcement
           pop(i).Schedule(Selected_behaviors) = pop(i).Schedule(Selected_behaviors) - pop(i).Schedule(Selected_behaviors) * (RR);
       end
       
       % Store Record for Current Iteration
       BestCosts(FEs) = BestSol.response; 
    
       %% ------- Rescheduling --------------------------------------------------    
       if std(pop(i).Schedule(1:N))== 0
           pop(i).Schedule = unifrnd(0.9,1,1,N);
           pop(i).Behaviors = LB+rand(1,N).*(UB-LB);
           %pop(i).response = feval('cec14_func',pop(i).Behaviors',ObjFunc_ID)  - (ObjFunc_ID*100);
           %pop(i).response = feval(fhd,pop(i).Behaviors',ObjFunc_ID); %CEC2017
           pop(i).response = Sphere(pop(i).Behaviors'); % Test Func
           disp(['-------------------------------- The  Learner ' num2str(i) ' is Rescheduled ']);
       end

   end % End for nPop
   
       
   %% Sort pop
   [~,SorteIndx] = sort([pop.response]);
   pop = pop(SorteIndx);
      
    
   % --- Show Iteration Information ---
   disp(['Iteration ' num2str(FEs) ': Best Cost = ' num2str(BestCosts(FEs))  ]);
  
      
 end % End While
BestSolCost=BestSol.response
end % Function PRO_v1()
function z=Sphere(x)

    z=sum(x.^2);
    
end