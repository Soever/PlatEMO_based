clc
clear
close all
%%
nPop=50; % 种群数

Max_iter=500; % 最大迭代次数

dim = 100; % 可选 2, 10, 30, 50, 100

%%  选择函数

Function_name=25; % 函数名： 1 - 30
[lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);

%% 调用算法
tic
[Best_score,Best_pos,cg_curve]=WOA(nPop,Max_iter,lb,ub,dim,fobj);
toc

%% plot
figure('Position',[400 200 300 250])
semilogy(cg_curve,'Color','r','Linewidth',1)
%     plot(cg_curve,'Color','r','Linewidth',1)
title(['Convergence curve, Dim=' num2str(dim)])
xlabel('Iteration');
ylabel(['Best score F' num2str(Function_name) ]);
axis tight
grid on
box on
set(gca,'color','none')
legend('WOA')

