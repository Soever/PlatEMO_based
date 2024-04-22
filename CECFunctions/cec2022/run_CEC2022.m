clc
clear
close all
%%
nPop=50; % 种群数
Max_iter=500; % 最大迭代次数
dim = 2; % 可选 2, 10, 20

%%  选择函数
Function_name=8; % 函数名： 1 - 12
[lb,ub,dim,fobj] = Get_Functions_cec2022(Function_name,dim);
% [Best_pos,Best_score,SO_curve]=SO(nPop,Max_iter,lb(1),ub(1),dim,fobj);
%% 调用算法
Optimal_results={}; % Optimal results
index = 1;

% WOA
tic
 

[Best_score,Best_pos,cg_curve]=WOA(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="WOA";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_pos;
Optimal_results{5,index}=toc;
index = index +1;
% HHO
tic
[Best_score,Best_pos,cg_curve]=WOA(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="HHO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_pos;
Optimal_results{5,index}=toc;

%% plot
figure
% semilogy(cg_curve,'Color','r','Linewidth',1)
for i = 1:size(Optimal_results, 2)
    plot(Optimal_results{2, i},'Linewidth',2)
    hold on
end
title(['Convergence curve, Dim=' num2str(dim)])
xlabel('Iteration');
ylabel(['Best score F' num2str(Function_name) ]);
axis tight
grid on
box on
set(gcf,'Position',[400 200 400 250])
legend(Optimal_results{1, :})

