clc ;
clear all ;

runs = 1 ;
T = 500  ;
N = 50 ;
dim = 10 ;

for i = 1: 10 
    Function_name=i; % 函数名： 1 - 12
    [lb,ub,dim,fobj] = Get_Functions_cec2022(Function_name,dim);
    lb = lb(1);
    ub = ub(1);
    for run = 1: runs 
        [Xfood, fval22,Convergence_curve,Trajectories,fitness_history, position_history]=SO_Cauchy(N,T,lb,ub,dim,fobj); 
        [Xfood, fval22,Convergence_curve,Trajectories,fitness_history, position_history]=SO_Gaussian(N,T,lb,ub,dim,fobj);
        [Xfood, fval22,Convergence_curve,Trajectories,fitness_history, position_history]=SO_T(N,T,lb,ub,dim,fobj);
    end
end
