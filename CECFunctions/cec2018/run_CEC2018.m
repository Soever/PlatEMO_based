% cec2018测试函数-关注微信公众号“优化算法侠”
clc;clear;close all
%%
nPop=50; % 种群数
Max_iter=300; % 最大迭代次数
dim = 10; % 维度，可选 2, 10, 20，30，50，100

%%  选择函数
for Function_name=1:30 % 函数名： 1 - 30

    % lb->下限，ub->上限，fobj->目标函数
    lb=-100*ones(1,dim);
    ub=100*ones(1,dim);
    fobj = @(x) cec18_func(x',Function_name); % 调用cec2018 测试函数
    
    %% 调用算法-（格式一样）
    Optimal_results={}; % 结果保存到 Optimal results
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
    % 如有更多算法，按上述格式在此处添加
    
    % 
    %% plot 收敛曲线
    figure(Function_name)
    for i = 1:size(Optimal_results, 2)
        %     plot(Optimal_results{2, i},'Linewidth',2)
        semilogy(Optimal_results{2, i},'Linewidth',2)
        hold on
    end
    title(['CEC2018 Convergence curve, Dim=' num2str(dim)])
    xlabel('Iteration');
    ylabel(['Best score on F' num2str(Function_name) ]);
    axis tight
    grid on;box on
    set(gcf,'Position',[400 200 400 250])
    legend(Optimal_results{1, :})
    saveas(gcf,['F' num2str(Function_name) '.jpg']) % 保存图窗为 jpg格式
    %% 保存到excel
    filename = ['cec2018_Results_D' num2str(dim) '.xlsx']; % 保存的文件名字
    sheet = 1; % 第一个sheet
    func =['F' num2str(Function_name)]; % 函数名字
    xlswrite(filename, [Optimal_results{1, :}], sheet, 'B1:D1') % 算法名字
    xlswrite(filename, {func}, sheet, ['A' num2str(Function_name+1)]) % 函数名字
    xlswrite(filename, [Optimal_results{3, :}], sheet, ['B' num2str(Function_name+1)]) % 最优结果

end