%% CEC2021 
clear all 
clc
SearchAgents_no=30; % Number of search agents 
%Éè¶¨²âÊÔº¯Êý
%Test functions are only defined for D=10, 9, 16, 18
%  F1 is defined on D=9
%  F2 is defined on D=16
%  F3 is defined on D=18
%  F4-F10 are defined on D=10.
Function_name=1; 
Max_iteration=500; % Maximum numbef of iterations 
D=10; 
lb=-100;
ub=100;

%cec21_basic_func;
% cec21_bias_func;
% cec21_shift_func;
% cec21_rot_func;
% cec21_bias_shift_func;
% cec21_bias_rot_func;
fobj = @(x) cec21_basic_func(x', Function_name);

[Best_pos,Best_score,SOA_curve]=SOA(SearchAgents_no,Max_iteration,lb,ub,D,fobj); %¿ªÊ¼ÓÅ»¯

%Draw objective space
figure
plot(SOA_curve,'Color','b','linewidth',1.5)
grid on;
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');
axis tight
grid on
box on
legend('SOA')
display(['The best solution obtained by SOA is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by SOA is : ', num2str(Best_score)]);