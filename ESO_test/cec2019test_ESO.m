clear all 
clc
runs = 30 ;
N=30;


% pool = parpool;
% 
resMean1 = zeros(10,1) ;
resStd1= zeros(10,1) ;
resMean2 = zeros(10,1) ;
resStd2= zeros(10,1) ;
parfor idx=1:10
    Function_idx=idx;
    [lb,ub,dim,fobj]= Get_Functions_cec2019(Function_idx);
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

 

