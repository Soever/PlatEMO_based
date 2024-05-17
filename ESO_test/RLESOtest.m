clear all 
clc
runs = 30 ;
N=30;

functionNames = arrayfun(@(x) ['F', num2str(x)], 1:23, 'UniformOutput', false);

T=500; 
stage = [1,2,3,4] ;%对应ESO四个改进的部分，数组中有哪些数字，则运行了哪些部分

Function_name=functionNames{20};
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
[Xfood1, fval1,Convergence_curve1,Trajectories1,fitness_history1, position_history1]=ESO_Algorithm(N,T,lb,ub,dim,fobj,stage);
[Xfood2, fval2,Convergence_curve,Trajectories2,fitness_history2, position_history2]=RLESO(N,T,lb,ub,dim,fobj,[1,3]); %[1,3]为动态参数部分





