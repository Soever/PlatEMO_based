clear all ;
clc

runs = 30 ;
T = 3000 ;
N = 100 ;
dim =30 ;%n==10||n==30||n==50||n==100)
function_num = 30;

fval1 = zeros(runs,function_num) ;
fval2 = zeros(runs,function_num) ;
fval3 = zeros(runs,function_num) ;
fval4 = zeros(runs,function_num) ;
result = zeros(function_num,4) ;
for i=1:function_num 
    Function_name=i; 

    [lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);
    lb = lb(1);
    ub = ub(1);
    
    disp(i);
    
   parfor run=1:runs
        [Xfood1, fval1(run,i),~]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
        [Xfood2, fval2(run,i),~]=MISO(N,T,lb,ub,dim,fobj); 
        % [Xfood2, fval2(run,i),~]=RLSO2_12(N,T,lb,ub,dim,fobj); 
        % [Xfood3, fval3(run,i),~]=RLSO2_13(N,T,lb,ub,dim,fobj); 
        %[Xfood4, fval4(run,i),~]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[1,2,3,4]);
    end
    
end
result(:,1)=mean(fval1,1);
result(:,2)=mean(fval2,1);
result(:,3)=mean(fval3,1);
result(:,4)=mean(fval4,1);
disp(result)
