clear all 
clc
runs = 30 ;
N=30;T=200; 


% pool = parpool;
% 
Function_name='F1';
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

results1= zeros(runs);
results2= zeros(runs);

[Xfood, fval,Convergence_curve,Trajectories,fitness_history, position_history]=ESO_1(N,T,lb,ub,dim,fobj); 

% parfor i=1:runs
%     [Xfood, fval,Convergence_curve,Trajectories,fitness_history, position_history]=ESO(N,T,lb,ub,dim,fobj); %¿ªÊ¼ÓÅ»¯
%     [Best_pos,Best_score,SO_curve]=SO(N,T,lb,ub,dim,fobj);
%     results1(i) = Best_score ;
%     results2(i) = fval ;
% end
% 
% resMean = mean(results1);
% resStd  = std(results1);
% 
% resMean2 = mean(results2);
% resStd2  = std(results2);