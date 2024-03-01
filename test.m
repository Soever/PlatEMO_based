clc
clear



N = 30 ;D = 3 ;
maxIteration = 500 ;runs = 30 ;n_func = 10 ;
% al_handles = {@SnakeOptimizer,@GA,@DE,@PSO};
% maxFE = maxIteration*N ; 
% for i = 1:n_func
%     func_name = "CEC2017_F"+ num2str(i); % 创建函数名字符串
%     func_handles{i} = str2func(func_name); % 创建函数句柄
% end
% obj_min= zeros(length(al_handles),n_func,runs);
% best = zeros(length(al_handles),n_func,runs);
% global result_t ;
% parfor i =1:length(al_handles)
%     for j = 1:n_func
%         for k = 1:runs
%             al =al_handles{i} ;
%             func = func_handles{j} ;
%             [~,obj,~] = platemo('problem',func,'algorithm',al,'N',N,'D',D,'maxFE',maxFE);
%             if i == 1 %SO通过全局变量result_t记录
%                 best(i,j,k) = result_t(end);
%             end
%             obj_min(i,j,k) = min(obj);
% 
%             % f = figure;
%             % plot(result_t,'LineWidth', 2, 'Color', 'k');
%             % title(strrep(char(al),'@','')+" on "+strrep(char(func),'@',''));
%         end
%     end
% end

rng(2024);

%platemo('problem',@CEC2017_F4,'algorithm',@GA,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2017_F1,'algorithm',@SnakeOptimizer,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2017_F1,'algorithm',@SnakeOptimizer_random,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2017_F1,'algorithm',@SnakeOptimizer_random2,'N',N,'D',D,'maxFE',15000);


% 
% platemo('problem',@CEC2017_F4,'algorithm',@PSO,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2017_F4,'algorithm',@SHADE,'N',N,'D',D,'maxFE',15000);

% runs = 5 ;
% Obj_min = zeros(runs,1) ;
% for i = 1:runs
%     [~,Obj,~] = platemo('problem',@CEC2017_F8,'algorithm',@GA,'N',N,'D',D,'maxFE',15000);
%     Obj_min = min(Obj);
% end
% display(mean(Obj_min));


