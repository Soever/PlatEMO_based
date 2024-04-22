clear all ;
clc

runs = 30 ;
T = 500 ;
N = 50 ;
dim =10 ;
function_num = 30 ;
result = zeros(function_num,8) ;
for i=1:function_num
    Function_name=i; 
    [lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);
    lb = lb(1);
    ub = ub(1);
    fval = zeros(8,runs) ;
    disp(i);
   parfor run=1:runs
        [Xfood1, fval1(run),Convergence_curve1]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
        [Xfood2, fval2(run),Convergence_curve1,Trajectories1,fitness_history1, position_history1]=SO_Cauchy(N,T,lb,ub,dim,fobj); 
        [Xfood2, fval3(run),Convergence_curve2,Trajectories2,fitness_history2, position_history2]=SO_Gaussian(N,T,lb,ub,dim,fobj);
        [Xfood3, fval4(run),Convergence_curve3,Trajectories3,fitness_history3, position_history3]=SO_T(N,T,lb,ub,dim,fobj);
        [Xfood2, fval5(run),Convergence_curve2,Trajectories2,fitness_history2, position_history2]=ESO_Gaussian(N,T,lb,ub,dim,fobj,[1,2,3,4]);
        [Xfood3, fval6(run),Convergence_curve3,Trajectories3,fitness_history3, position_history3]=ESO_T(N,T,lb,ub,dim,fobj,[1,2,3,4]);
        [Xfood3, fval7(run),Convergence_curve3,Trajectories3,fitness_history3, position_history3]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[2]);
        [Xfood3, fval8(run),Convergence_curve3,Trajectories3,fitness_history3, position_history3]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[1,2,3,4]);
    end
    result(i,1)=mean(fval1);
    result(i,2)=mean(fval2);
    result(i,3)=mean(fval3);
    result(i,4)=mean(fval4);
    result(i,5)=mean(fval5);
    result(i,6)=mean(fval6);
    result(i,7)=mean(fval7);
    result(i,8)=mean(fval8);
end

disp(result)
