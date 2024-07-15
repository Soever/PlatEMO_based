clear all ;
clc


runs =51 ;
T =2000;

N = 100 ;
dim =20 ;
function_num = 10 ;

fval1 = zeros(runs,function_num) ;
fval2 = zeros(runs,function_num) ;
fval3 = zeros(runs,function_num) ;
fval4 = zeros(runs,function_num) ;

fval5 = zeros(runs,function_num) ;
result = zeros(function_num,5) ;
for i=5:5

    
Function_name = i ;


    [lb,ub,dim,fobj] = Get_Functions_cec2019(Function_name);
    T = dim*10000/N;
    lb = lb(1);
    ub = ub(1);
    
    disp(i);
    
   parfor run=1:runs
        [Xfood1, fval1(run,i),~]=SO_Algorithm(N,T,lb,ub,dim,fobj); 

        %[Xfood1, fval5(run,i),~]=RLSO2_14(N,T,lb,ub,dim,fobj); 
        [Xfood2, fval2(run,i),~]=RLSO5_1(N,T,lb,ub,dim,fobj); 
        %[Xfood3, fval3(run,i),~]=RLSO2_18(N,T,lb,ub,dim,fobj); 
        [Xfood3, fval3(run,i),~]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[1,2,3,4]);
        [Xfood4, fval4(run,i),~]=ESO_Algorithm(N,T/2,lb,ub,dim,fobj,[1,2,3,4]);
    end

    
end
result(:,1)=mean(fval1,1);
result(:,2)=mean(fval2,1);
result(:,3)=mean(fval3,1);
result(:,4)=mean(fval4,1);
result(:,5)=mean(fval5,1);
disp(result)
