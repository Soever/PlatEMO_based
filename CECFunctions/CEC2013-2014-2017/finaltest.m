
clc
clear all ;
runs =1;
T =30; % 最大迭代次数
N =100; % 种群数
dim =10; % 维度

rng(2024);

function_num = 30 ;
fval1 = zeros(runs,function_num) ;
fval2 = zeros(runs,function_num) ;
fval3 = zeros(runs,function_num) ;
fval4 = zeros(runs,function_num) ;
fval5 = zeros(runs,function_num) ;
fval6 = zeros(runs,function_num) ;
fval7 = zeros(runs,function_num) ;
fval8 = zeros(runs,function_num) ;
fval9 = zeros(runs,function_num) ;
fval10 = zeros(runs,function_num) ;
fval11 = zeros(runs,function_num) ;
fval12 = zeros(runs,function_num) ;
fval13 = zeros(runs,function_num) ;
fval14 = zeros(runs,function_num) ;
result2014 = zeros(function_num,14) ;
curve1 = zeros(runs, function_num, T);
curve2 = zeros(runs, function_num, T);
curve3 = zeros(runs, function_num, T);
curve4 = zeros(runs, function_num, T);
curve5 = zeros(runs, function_num, T);
curve6 = zeros(runs, function_num, T);
curve7 = zeros(runs, function_num, T);
curve8 = zeros(runs, function_num, T);
curve9 = zeros(runs, function_num, T);
curve10 = zeros(runs, function_num, T);
curve11 = zeros(runs, function_num, T);
curve12 = zeros(runs, function_num, T);

for i=1:30
    Function_name=i; 

    [lb,ub,dim,fobj] = Get_Functions_cec2014(Function_name,dim);
    % lb = lb(1);
    % ub = ub(1);

    disp(i);

   for run=1:runs
       [fval1(run,i), Best_pos1, curve1(run,i,:)] = GWO(N, T, lb, ub, dim, fobj);
[fval2(run,i), Best_pos2, curve2(run,i,:)] = HHO(N, T, lb, ub, dim, fobj);
[fval3(run,i), Best_pos3, curve3(run,i,:)] = SCA(N, T, lb, ub, dim, fobj);
[fval4(run,i), Best_pos4, curve4(run,i,:)] = GJO(N, T, lb, ub, dim, fobj);
[Best_pos5, fval5(run,i), curve5(run,i,:)] = SOA(N, T, lb, ub, dim, fobj);
[Best_pos6,fval6(run,i),  curve6(run,i,:)] = SO_Algorithm(N, T, lb(1), ub(1), dim, fobj);
[fval7(run,i), Best_pos7, curve7(run,i,:)] = AGWO(N, T, lb(1), ub(1), dim, fobj);
[Best_pos8,fval8(run,i),  curve8(run,i,:),~,~,~] = ESO_AlgorithmFE(N, T/2, lb(1), ub(1), dim, fobj, [1, 2, 3, 4]);
[Best_pos9,fval9(run,i),  curve9(run,i,:)] = fdb_tlabc_me(N, T, dim, lb, ub, fobj);
[fval10(run,i), Best_pos10, curve10(run,i,:)] = PDWOA_me(N, T, lb, ub, dim, fobj);
[fval11(run,i), Best_pos11, curve11(run,i,:), ~] = QQLMPA_me(N, T, lb, ub, dim, fobj);
[Best_pos12,fval12(run,i),  curve12(run,i,:)] = RLSO5_3(N, T, lb(1), ub(1), dim, fobj);

      % [Xfood1, fval1(run,i),~]=RLSO5_6(N,T,lb,ub,dim,fobj); 
      %  [Xfood2, fval2(run,i),~]=ESO_AlgorithmFE(N,T/2,lb,ub,dim,fobj,[1,2,3,4]); 
        %[Xfood3, fval3(run,i),~,~,~, ~]=RLSO2_1(N,T,lb,ub,dim,fobj); 
       % [Xfood3, fval3(run,i),~]=RLSO5_4(N,T,lb,ub,dim,fobj); 
       % [Xfood4, fval4(run,i),~]=RLSO5_5(N,T,lb,ub,dim,fobj);
        
        %[Xfood4, fval4(run,i),~,~,~, ~]=ESO_Algorithm(N,T,lb,ub,dim,fobj,[1,2,3,4]);
    end

end
result2014(:, 1) = mean(fval1, 1);
result2014(:, 2) = mean(fval2, 1);
result2014(:, 3) = mean(fval3, 1);
result2014(:, 4) = mean(fval4, 1);
result2014(:, 5) = mean(fval5, 1);
result2014(:, 6) = mean(fval6, 1);
result2014(:, 7) = mean(fval7, 1);
result2014(:, 8) = mean(fval8, 1);
result2014(:, 9) = mean(fval9, 1);
result2014(:, 10) = mean(fval10, 1);
result2014(:, 11) = mean(fval11, 1);
result2014(:, 12) = mean(fval12, 1);
