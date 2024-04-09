% 微信公众号：优化算法侠，Swarm-Opti
clc;clear;close all
%%
nPop=50; % 种群数
Max_iter=500; % 最大迭代次数

%%  选择函数
dim = 20; % 维数，可选 2, 5, 10, 15, 20
Function_name=10; % 函数名： 1 - 10
% lb->下限，ub->上限，fobj->目标函数，dim-> 维度
[lb,ub,dim,fobj] = Get_Functions_cec2020(Function_name, dim);

%% 调用算法
Optimal_results={}; % 保存Optimal results
index = 1;
% PSO
tic
[Best_score,Best_x,cg_curve]=PSO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="PSO";           % 算法名字
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
% GTO
tic
[Best_score,Best_x,cg_curve]=GTO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="GTO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
% AVOA
tic
[Best_score,Best_x,cg_curve]=AVOA(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="AVOA";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
% SHO
tic
[Best_score,Best_x,cg_curve]=SHO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="SHO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;

%% plot 绘图
figure

for i = 1:size(Optimal_results, 2)
%     plot(Optimal_results{2, i},'Linewidth',2)
    semilogy(Optimal_results{2, i},'Linewidth',2)
    hold on
end
title(['Convergence curve, Dim=' num2str(dim)])
xlabel('Iteration');
ylabel(['Best score F' num2str(Function_name) ]);
ax = gca;
set(ax,'Tag',char([100,105,115,112,40,39,20316,32773,58,...
    83,119,97,114,109,45,79,112,116,105,39,41]));
axis tight
grid on
box on
eval(ax.Tag)
set(gcf,'Position',[400 200 400 250])
legend(Optimal_results{1, :})

