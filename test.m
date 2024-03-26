clc
clear
N = 100 ;D = 2 ;
maxIteration = 1000 ;runs = 30 ;n_func = 10 ;
maxFE = maxIteration*N ;

% func_handles = get_benchmark_func("ALL");
% 
% Tc = {0.4,0.5,0.6,0.7,0.8,0.9};
% Pc = {1.2,1.5,1.8,2,2.3,2.6} ;
% al_handles = {} ;
% for i = 1:length(Tc)
%     for j = i:length(Pc)
%         al_handles{end+1} = {@SnakeOptimizer,Pc{j},Tc{i},0.05,2};
%     end
% end


global result_t ;
% al_handles={
%     {@SnakeOptimizer,0.15,0.6,0.5,0.05,2};
%     {@SnakeOptimizer,0.2,0.6,0.5,0.05,2};
%     {@SnakeOptimizer,0.25,0.6,0.5,0.05,2};
%     {@SnakeOptimizer,0.3,0.6,0.5,0.05,2};
%     {@SnakeOptimizer,0.35,0.6,0.5,0.05,2};
% };



% obj_min= zeros(length(al_handles),n_func,runs);
% best = zeros(length(al_handles),n_func,runs);


% for i =1:length(al_handles)
%     parfor j = 1:length(func_handles)
%         for k = 1:10
%             al =al_handles{i} ;
%             func = func_handles{j} ;
%             [~,obj,~] = platemo('problem',func,'algorithm',al,'N',N,'maxFE',maxFE);
%             if i == 1 %SO通过全局变量result_t记录
%                 best(i,j,k) = result_t(end);
%             end
%             obj_min(i,j,k) = min(obj);
%         end
%     end
% end



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
platemo('problem',@CEC2008_F1,'algorithm',@ESO1 ,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2008_F1,'algorithm',@ESO2 ,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2008_F1,'algorithm',@ESO3 ,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2008_F1,'algorithm',@ESO4 ,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2008_F1,'algorithm',@SnakeOptimizer ,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2008_F1,'algorithm',@eso_emo ,'N',N,'D',D,'maxFE',15000);

% platemo('problem',@CEC2008_F1,'algorithm',@SO1,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2008_F1,'algorithm',@SO2,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2008_F1,'algorithm',@SO3,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2008_F1,'algorithm',@SO13,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2008_F1,'algorithm',@SO23,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2008_F1,'algorithm',@SO12,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2020_F1,'algorithm',@LA_PSO,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2020_F1,'algorithm',@PSO,'N',N,'D',D,'maxFE',15000);
%platemo('problem',@CEC2017_F4,'algorithm',@GA,'N',N,'D',D,'maxFE',15000);
%platemo('problem',@CEC2017_F1,'algorithm',@SnakeOptimizer,'N',N,'D',D,'maxFE',15000);
%platemo('problem',@CEC2017_F1,'algorithm',@SnakeOptimizer_pro,'N',N,'D',D,'maxFE',15000);
%platemo('problem',@ZDT1,'algorithm',@MOPSO,'N',N,'D',D,'maxFE',15000);
%platemo('problem',@ZDT1,'algorithm',@MRPSO,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2017_F1,'algorithm',@SnakeOptimizer_random,'N',N,'D',D,'maxFE',15000);
% platemo('problem',@CEC2017_F1,'algorithm',@SnakeOptimizer_random2,'N',N,'D',D,'maxFE',15000);


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


