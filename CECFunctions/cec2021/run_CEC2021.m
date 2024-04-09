clc
clear
close all
%%
nPop=50; % ��Ⱥ��
Max_iter=500; % ����������
dim = 20; % ά�ȣ���ѡ 2, 10, 20

%%  ѡ����
Function_name=10; % �������� 1 - 10
% lb->���ޣ�ub->���ޣ�fobj->Ŀ�꺯��
[lb,ub,dim,fobj] = Get_Functions_cec2021(Function_name,dim);

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
