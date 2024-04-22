clear all 
clc
runs = 30 ;
N=100;



Function_name='F1';
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

results1= zeros(runs);
results2= zeros(runs);
 a = [1,2,3,4];
rng(2024);
 T=1000; 
%[Xfood, fval22,Convergence_curve]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
% [Xfood, fval22,Convergence_curve]=SO_Algorithm(N,T,lb,ub,dim,fobj); 
%[Best_pos,Best_score,SO_curve]=RLSO_FitnessDiversityGrid(N,T,lb,ub,dim,fobj); 
  [Xfood, fval22,Convergence_curve,Trajectories,fitness_history, position_history]=RLSO2_4(N,T,lb,ub,dim,fobj); 
  % [Xfood, fval22,Convergence_curve,Trajectories,fitness_history, position_history]=ESO(N,T,lb,ub,dim,fobj); 
% [Best_pos,Best_score,SO_curve]=RLSO_FitnessDiversityGrid(N,T,lb,ub,dim,fobj); 
  [Best_pos,Best_score,SO_curve]=SO_Algorithm(N,T,lb,ub,dim,fobj);
 % [Best_pos,Best_score,SO_curve,~,~]=RLSO_2(N,T,lb,ub,dim,fobj);
 figure('Position',[39         479        1727         267])
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
title('Convergence curve')
xlabel('Iteration#');
ylabel('Best score obtained so far');
box on
legend('RLSO_2','SO')




