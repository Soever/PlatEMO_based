clear all 
clc
runs = 30 ;
N=100;



Function_name='F20';
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

results1= zeros(runs);
results2= zeros(runs);
 a = [1,2,3,4]; 
% rng(2024);
 T=1000; 
%[Xfood, fval22,Convergence_curve]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
% [Xfood, fval22,Convergence_curve]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
%[Best_pos,Best_score,SO_curve]=RLSO_FitnessDiversityGrid(N,T,lb,ub,dim,fobj); 
   %[Xfood, fval22,Convergence_curve,Trajectories,fitness_history, position_history]=RLSO2_11(N,T,lb,ub,dim,fobj);
   [Xfood, fval22,Convergence_curve]=RLSO2_14(N,T,lb,ub,dim,fobj);
     % [Xfood, fval22,Convergence_curve,~,~, ~,fitness_history1]=SO_operator_compare(N,T,lb,ub,dim,fobj,1);
     % [Xfood, fval22,Convergence_curve,~,~, ~,fitness_history2]=SO_operator_compare(N,T,lb,ub,dim,fobj,2);
     % [Xfood, fval22,Convergence_curve,~,~, ~,fitness_history3]=SO_operator_compare(N,T,lb,ub,dim,fobj,3);
     % [Xfood, fval22,Convergence_curve,~,~, ~,fitness_history4]=SO_operator_compare(N,T,lb,ub,dim,fobj,4);
     % [Xfood, fval22,Convergence_curve,fitness_history]=so_oo(N,T,lb,ub,dim,fobj);
    % [Best_pos,Best_score,SO_curve]=so_3(N,T,lb,ub,dim,fobj);

% hold on
% title(Function_name);
% semilogy(min(fitness_history1,[],1),'Color','m','linewidth',1.5)
% semilogy(min(fitness_history2,[],1),'Color','y','linewidth',1.5)
% semilogy(min(fitness_history3,[],1),'Color','b','linewidth',1.5)
% semilogy(min(fitness_history4,[],1),'Color','g','linewidth',1.5)
% semilogy(min(fitness_history,[],1),'Color','r','linewidth',1.5)
% legend("1","2","3","4","SO");
% plot(min(fitness_history1,[],1)) ;
% [Best_pos,Best_score,SO_curve]=RLSO2_random(N,T,lb,ub,dim,fobj);
    % [Xfood, fval22,Convergence_curve,Trajectories,fitness_history, position_history]=ESO(N,T,lb,ub,dim,fobj); 
% [Best_pos,Best_score,SO_curve]=RLSO_FitnessDiversityGrid(N,T,lb,ub,dim,fobj); 
   %%
   %[Best_pos,Best_score,SO_curve]=RLSO4_1(N,T,lb,ub,dim,fobj);
   %[Best_pos,Best_score,SO_curve]=SO_Algorithm(N,T,lb,ub,dim,fobj);
   [Best_pos1,Best_score1,ESO_curve]=ESO_AlgorithmFE(N,T/2,lb,ub,dim,fobj,[1,2,3,4]);
  % [Best_pos1,Best_score1,ESO_curve,~,~]=RLSO2_15(N,T,lb,ub,dim,fobj);
   [Best_pos1,Best_score1,SO_curve]=RLSO2_16(N,T,lb,ub,dim,fobj);
 % [Best_pos,Best_score,SO_curve,~,~]=RLSO_2(N,T,lb,ub,dim,fobj);
 fig = figure;

% 获取屏幕的大小
screenSize = get(0, 'ScreenSize'); % 返回形如[left, bottom, width, height]
screenWidth = screenSize(3);
screenHeight = screenSize(4);

% 设定图形窗口的宽度和高度，保持1:2的比例
figWidth = screenWidth / 1.2;  % 例如，设置图形窗口宽度为屏幕宽度的1/4
figHeight =  figWidth/2.5;    % 高度是宽度的两倍

% 计算图形窗口的左下角坐标，以使图形窗口居中
positionLeft = (screenWidth - figWidth) / 2;
positionBottom = (screenHeight - figHeight) / 2;

% 设置图形窗口的位置和大小
set(fig, 'Position', [positionLeft, positionBottom, figWidth, figHeight]);

color1 = [205 205 0];
color2 = [139 101 8];
color3 = [205 155 155];
color4 = [238 121 66];
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
box on
axis tight

subplot(1,2,2);
semilogy(Convergence_curve,'Color','r','linewidth',1.5)
hold on
semilogy(SO_curve,'Color','b','linewidth',1.5)
semilogy(ESO_curve,'Color','y','linewidth',1.5)
title('Convergence curve')
xlabel('Iteration#');
ylabel('Best score obtained so far');
box on
%legend('RLSO_{2}','RLSO_{random}')
legend('RLSO_{2}','SO','ESO');
%legend('SO_{234}','SO');


