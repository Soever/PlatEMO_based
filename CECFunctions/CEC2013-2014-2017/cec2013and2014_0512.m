
clc
clear all ;
runs =1;
T =30;
N =100;
dim =50;
% function_num = 28 ;

% fval1 = zeros(runs,function_num) ;
% fval2 = zeros(runs,function_num) ;
% fval3 = zeros(runs,function_num) ;
% fval4 = zeros(runs,function_num) ;
% result2013 = zeros(function_num,4) ;
% 
% 
% for i=1:function_num
%     Function_name=i; 
% 
%     [lb,ub,dim,fobj] = Get_Functions_cec2013(Function_name,dim);
%     lb = lb(1);
%     ub = ub(1);
% 
%     disp(i);
% 
%    parfor run=1:runs
%         [Xfood1, fval1(run,i),~]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
%         [Xfood2, fval2(run,i),~,~,~, ~]=RLSO2_11(N,T,lb,ub,dim,fobj); 
%         [Xfood3, fval3(run,i),~,~,~, ~]=RLSO2_1(N,T,lb,ub,dim,fobj); 
%         [Xfood4, fval4(run,i),~,~,~, ~]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[1,2,3,4]);
%     end
% 
% end
% result2013(:,1)=mean(fval1,1);
% result2013(:,2)=mean(fval2,1);
% result2013(:,3)=mean(fval3,1);
% result2013(:,4)=mean(fval4,1);
% 
% disp(result2013)

function_num = 30 ;
fval1 = zeros(runs,function_num) ;
fval2 = zeros(runs,function_num) ;
fval3 = zeros(runs,function_num) ;
fval4 = zeros(runs,function_num) ;
fval5 = zeros(runs,function_num) ;
result2014 = zeros(function_num,5) ;

for i=29:30
    Function_name=i; 

    [lb,ub,dim,fobj] = Get_Functions_cec2014(Function_name,dim);
    lb = lb(1);
    ub = ub(1);

    disp(i);

   for run=1:runs
       % [Xfood1, fval1(run,i),~]=RLSO5_6(N,T,lb,ub,dim,fobj); 
      %  [Xfood2, fval2(run,i),~]=ESO_AlgorithmFE(N,T/2,lb,ub,dim,fobj,[1,2,3,4]); 
        %[Xfood3, fval3(run,i),~,~,~, ~]=RLSO2_1(N,T,lb,ub,dim,fobj); 
       % [Xfood3, fval3(run,i),~]=RLSO5_4(N,T,lb,ub,dim,fobj); 
       % [Xfood4, fval4(run,i),~]=RLSO5_5(N,T,lb,ub,dim,fobj);
        [Xfood5, fval5(run,i),~]=RLSO5_7(N,T,lb,ub,dim,fobj);
        %[Xfood4, fval4(run,i),~,~,~, ~]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[1,2,3,4]);
    end

end
result2014(:,1)=mean(fval1,1);
result2014(:,2)=mean(fval2,1);
result2014(:,3)=mean(fval3,1);
result2014(:,4)=mean(fval4,1);
result2014(:,5)=mean(fval5,1);
