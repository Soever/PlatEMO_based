clc
clear
close all
%%
nPop=50; % 种群数
Max_iter=500; % 最大迭代次数
dim = 20; % 维度，可选 2, 10, 20

%%  选择函数
Function_name=10; % 函数名： 1 - 10
% lb->下限，ub->上限，fobj->目标函数
[lb,ub,dim,fobj] = Get_Functions_cec2021(Function_name,dim);

%% 调用算法
Optimal_results={}; % 保存Optimal results
index = 1;
% WOA
tic
[Best_score,Best_x,cg_curve]=WOA(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="WOA";         % 算法名字
Optimal_results{2,index}=cg_curve;      % 收敛曲线
Optimal_results{3,index}=Best_score;   % 最优函数值
Optimal_results{4,index}=Best_x;          % 最优变量
Optimal_results{5,index}=toc;               % 运行时间
index = index +1;
% HHO
tic
[Best_score,Best_x,cg_curve]=HHO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="HHO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
% GWO
tic
[Best_score,Best_x,cg_curve]=GWO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="GWO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
%% plot
figure

for i = 1:size(Optimal_results, 2)
%     plot(Optimal_results{2, i},'Linewidth',2)
    semilogy(Optimal_results{2, i},'Linewidth',2)
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

