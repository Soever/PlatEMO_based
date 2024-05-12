clc
clear
close all
%%
nPop=50; % 种群数

Max_iter=500; % 最大迭代次数

dim = 20; % 可选
runs = 50 ;
%%  选择函数

for Function_name = 1:12
    [lb,ub,dim,fobj] = Get_Functions_cec2017(Function_name,dim);
    lb = lb(1);ub = ub(1);
    for run = 1:runs
        [Best_score1,Best_pos1,cg_curve1,~,~]=RLSO2_1(nPop,Max_iter,lb,ub,dim,fobj);
        [Best_score1,Best_pos1,cg_curve1,~,~]=RLSO2_11(nPop,Max_iter,lb,ub,dim,fobj);
        [Best_score2,Best_pos2,cg_curve2]=SO_Algorithm(nPop,Max_iter,lb,ub,dim,fobj);
    end
end
Function_name=25; % 函数名： 1 - 30


%% 调用算法
tic

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

