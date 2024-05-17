function [outputArg1,outputArg2] = test24(inputArg1,inputArg2)
clear all 
clc
runs = 30 ;
N=30;

functionNames = arrayfun(@(x) ['F', num2str(x)], 1:23, 'UniformOutput', false);

resMean =zeros(23,6);
resStd =zeros(23,6);


T=50; 
parfor index=1:23
    Function_name=functionNames{index};
    [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
    results1= zeros(runs,1);
    results2= zeros(runs,1);
    results3= zeros(runs,1);
    results4= zeros(runs,1);
    resultsSO= zeros(runs,1);
    resultsESO= zeros(runs,1);
    for i=1:runs
        [Xfood, fval,Convergence_curve,Trajectories,fitness_history, position_history]=ESO(N,T,lb,ub,dim,fobj); 
       
        [Xfood1, fval1,Convergence_curve1,Trajectories1,fitness_history1, position_history1]=ESO_1(N,T,lb,ub,dim,fobj); 
        [Xfood2, fval2,Convergence_curve2,Trajectories2,fitness_history2, position_history2]=ESO_2(N,T,lb,ub,dim,fobj); 
        [Xfood3, fval3,Convergence_curve3,Trajectories3,fitness_history3, position_history3]=ESO_3(N,T,lb,ub,dim,fobj); 
        [Xfood4, fval4,Convergence_curve4,Trajectories4,fitness_history4, position_history4]=ESO_4(N,T,lb,ub,dim,fobj); 
        [Best_pos,Best_score,SO_curve]=SO(N,T,lb,ub,dim,fobj);
        resultsESO(i) = fval ;
        results1(i) = fval1 ;
        results2(i) = fval2 ;
        results3(i) = fval3 ;
        results4(i) = fval4 ;
        resultsSO(i) = Best_score ;
        
    end
     resMean(index, :) = [mean(resultsSO), mean(results1), mean(results2), mean(results3), mean(results4), mean(resultsESO)];
     resStd(index, :) = [std(resultsSO), std(results1), std(results2), std(results3), std(results4), std(resultsESO)];

    % resMean(index,1) = mean(resultsSO);
    % resStd (index,1) = std(resultsSO);
    % resMean(index,2) = mean(results1);
    % resStd (index,2) = std(results1);
    % resMean(index,3) = mean(results2);
    % resStd (index,3) = std(results2);
    % resMean(index,4) = mean(results3);
    % resStd (index,4) = std(results3);
    % resMean(index,5) = mean(results4);
    % resStd (index,5) = std(results4);
    % resMean(index,6) = mean(resultsESO);
    % resStd (index,6) = std(resultsESO);
end
%%
clear all 
clc
runs = 30 ;
N=30;T=200; 

Function_name='cec06';
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);



results1= zeros(runs,1);
results2= zeros(runs,1);
parfor i=1:runs
    [Xfood, fval,Convergence_curve,Trajectories,fitness_history, position_history]=ESO(N,T,lb,ub,30,fobj); %¿ªÊ¼ÓÅ»¯
    [Best_pos,Best_score,SO_curve]=SO(N,T,lb,ub,30,fobj);
    results1(i) = Best_score ;
    results2(i) = fval ;
end

resMean = mean(results1);
resStd  = std(results1);

resMean2 = mean(results2);
resStd2  = std(results2);



a = 5 ;
end

