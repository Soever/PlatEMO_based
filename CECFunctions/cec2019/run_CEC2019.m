% ΢�Ź��ںţ��Ż��㷨����Swarm-Opti
clc;clear;close all
%%
nPop=50; % ��Ⱥ��
Max_iter=500; % ����������

%%  ѡ����
Function_name=10; % �������� 1 - 10
% lb->���ޣ�ub->���ޣ�fobj->Ŀ�꺯����dim-> ά��
[lb,ub,dim,fobj] = Get_Functions_cec2019(Function_name);

%% �����㷨
Optimal_results={}; % ����Optimal results
index = 1;
% WOA
tic
[Best_score,Best_x,cg_curve]=WOA(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="WOA";         % �㷨����
Optimal_results{2,index}=cg_curve;      % ��������
Optimal_results{3,index}=Best_score;   % ���ź���ֵ
Optimal_results{4,index}=Best_x;          % ���ű���
Optimal_results{5,index}=toc;               % ����ʱ��
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
% DBO
tic
[Best_score,Best_x,cg_curve]=DBO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="DBO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
% COA
tic
[Best_score,Best_x,cg_curve]=COA(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="COA";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
% PSO
tic
[Best_score,Best_x,cg_curve]=PSO(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="PSO";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
% YDSE
tic
[Best_score,Best_x,cg_curve]=YDSE(nPop,Max_iter,lb,ub,dim,fobj);
Optimal_results{1,index}="YDSE";
Optimal_results{2,index}=cg_curve;
Optimal_results{3,index}=Best_score;
Optimal_results{4,index}=Best_x;
Optimal_results{5,index}=toc;
index = index +1;
%% plot ��ͼ
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

