clear all ;
clc

runs = 30 ;
T = 500 ;
N = 100 ;
dim =100 ;
function_num = 30 ;
result = zeros(function_num,2) ;
for i=1:function_num
    Function_name=i; 
    [lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);
    lb = lb(1);
    ub = ub(1);
    fval = zeros(8,runs) ;
    disp(i);
    
   parfor run=1:runs
        [Xfood1, fval1(run),Convergence_curve1]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
        [Xfood2, fval2(run),Convergence_curve1,Trajectories1,fitness_history1, position_history1]=RLSO2_11(N,T,lb,ub,dim,fobj); 
    end
    result(i,1)=mean(fval1);
    result(i,2)=mean(fval2);

end

disp(result)
