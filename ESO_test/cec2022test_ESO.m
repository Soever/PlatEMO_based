clear all 
clc
runs = 30 ;
N=30;


% pool = parpool;
% 
fun_num = 12 ;
resMean1 = zeros(fun_num,1) ;
resStd1= zeros(fun_num,1) ;
resMean2 = zeros(fun_num,1) ;
resStd2= zeros(fun_num,1) ;
parfor idx=1:fun_num
    Function_idx=idx;
    [lb,ub,dim,fobj]= Get_Functions_cec2022(Function_idx,20);
    lb = lb(1);
    ub = ub(1);
    results1= zeros(runs,1);
    results2= zeros(runs,1);
    a = [1,2,3,4];
    
    T=1000; 
    for i=1:runs
        [Xfood, fval,Convergence_curve,Trajectories,fitness_history, position_history]=ESO_Algorithm(N,T,lb,ub,dim,fobj,a); 
        [Best_pos,Best_score,SO_curve]=SO(N,T,lb,ub,dim,fobj); 
        results1(i) =fval  ;
         results2(i) =Best_score;
    end
    resMean1(idx) = mean(results1);
    resMean2(idx) = mean(results2);
    resStd1(idx) = std(results1);
    resStd2(idx) = std(results2);
    % 
end

 

