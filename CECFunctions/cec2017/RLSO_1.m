clear all ;
clc

runs = 100 ;
T = 500 ;
N = 50 ;
dim =10 ;
function_num = 30 ;
result = zeros(function_num,1) ;
for i=1:function_num
    Function_name=i; 
    [lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);
    lb = lb(1);
    ub = ub(1);
    fval = zeros(1,runs) ;
   parfor run=1:runs
        [Xfood1, fval1(run),Convergence_curve1]=RLSO_FitnessDiversityGrid(N,T,lb,ub,dim,fobj); 
   end
    result(i,1)=mean(fval1);
end

disp(result)
