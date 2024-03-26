
clear all 
clc

N=30; % Number of search agents

Function_name='F18'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

T=200; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
  [Xfood, fval,Convergence_curve,Trajectories,fitness_history, position_history]=ESO(N,T,lb,ub,dim,fobj); %¿ªÊ¼ÓÅ»¯
  [Best_pos,Best_score,SO_curve]=SO(N,T,lb,ub,dim,fobj); 
 figure('Position',[39         479        1727         267])
color1 = [205 205 0];
color2 = [139 101 8];
color3 = [205 155 155];
color4 = [238 121 66];
%Draw search space
subplot(1,5,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])
box on
axis tight

subplot(1,5,2);
semilogy(Convergence_curve,'Color','r','linewidth',1.5)
hold on
semilogy(SO_curve,'Color','b','linewidth',1.5)
title('Convergence curve')
xlabel('Iteration#');
ylabel('Best score obtained so far');
box on
legend('ESO','SO')

axis tight

subplot(1,5,3);
hold on
semilogy(Trajectories(1,:),'Color',color4/255,'linewidth',1.5);
title('Trajectory ')
xlabel('Iteration#')
box on
axis tight

subplot(1,5,4);
hold on
a=mean(fitness_history);
semilogy(a,'Color',color2/255,'linewidth',1.5);
title('Average Fitness ')
xlabel('Iteration#')
box on
axis tight


subplot(1,5,5);
hold on
for k1 = 1: size(position_history,1)
    for k2 = 1: size(position_history,2)
        plot(position_history(k1,k2,1),position_history(k1,k2,2),'.','markersize',1,'MarkerEdgeColor','k','markerfacecolor','k');
    end
end
plot(Xfood(1),Xfood(2),'.','markersize',10,'MarkerEdgeColor','r','markerfacecolor','r','linewidth',2);
title('Search history (x1 and x2 only)')
xlabel('x1')
ylabel('x2')
box on
axis tight
% 
subplot(1,5,5);
hold on
func_plot1(Function_name)
