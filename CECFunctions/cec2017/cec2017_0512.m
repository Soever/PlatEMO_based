clear all ;
clc

runs =30;

N = 100 ;
dim =100 ;%n==10||n==30||n==50||n==100)
T = dim*10000/N ;
function_num = 30;

fval1 = zeros(runs,function_num) ;
fval2 = zeros(runs,function_num) ;
fval3 = zeros(runs,function_num) ;
fval4 = zeros(runs,function_num) ;
fval5 = zeros(runs,function_num) ;
fval6 = zeros(runs,function_num) ;
fval7 = zeros(runs,function_num) ;
result = zeros(function_num,7) ;
for i=1:function_num 
    Function_name=i; 

    [lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);
    lb = lb(1);
    ub = ub(1);
    
    disp(i);
    
   parfor run=1:runs
        [Xfood1, fval1(run,i),~]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
        [Xfood2, fval2(run,i),~,~,~]=ESO_Algorithm(N,T/2,lb,ub,dim,fobj,[1,2,3,4]); 
        [Xfood3, fval3(run,i),~]=RLSO5_3(N,T,lb,ub,dim,fobj); 
        [Xfood4, fval4(run,i),~]=RLSO5_4(N,T,lb,ub,dim,fobj); 
        [Xfood5, fval5(run,i),~]=RLSO5_5(N,T,lb,ub,dim,fobj); 
        [Xfood6, fval6(run,i),~]=RLSO5_6(N,T,lb,ub,dim,fobj); 
        [Xfood7, fval7(run,i),~]=RLSO5_7(N,T,lb,ub,dim,fobj);  
        % [Xfood2, fval2(run,i),~]=RLSO2_12(N,T,lb,ub,dim,fobj); 
        % [Xfood3, fval3(run,i),~]=RLSO2_13(N,T,lb,ub,dim,fobj); 
        %[Xfood4, fval4(run,i),~]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[1,2,3,4]);
    end
    
end
result(:,1)=mean(fval1,1);
result(:,2)=mean(fval2,1);
result(:,3)=mean(fval3,1);
result(:,4)=mean(fval4,1);
result(:,5)=mean(fval5,1);
result(:,6)=mean(fval6,1);
result(:,7)=mean(fval7,1);
disp(result)
