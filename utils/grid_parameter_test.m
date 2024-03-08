clc
clear
N = 100 ;D = 2 ;
maxIteration = 1000 ;runs = 30 ;n_func = 10 ;
maxFE = maxIteration*N ;

func_handles = get_benchmark_func("CEC2008");

Tc = {0.4,0.5,0.6,0.7,0.8,0.9};
Pc = {1.2,1.5,1.8,2,2.3,2.6} ;
al_handles = {} ;
for i = 1:length(Tc)
    for j = i:length(Pc)
        al_handles{end+1} = {@SnakeOptimizer,Pc{j},Tc{i},0.05,2};
    end
end


global result_t ;


obj_min= zeros(length(al_handles),n_func,runs);
best = zeros(length(al_handles),n_func,runs);


parfor i =1:length(al_handles)
    for j = 1:length(func_handles)
        for k = 1:10
            al =al_handles{i} ;
            func = func_handles{j} ;
            [~,obj,~] = platemo('problem',func,'algorithm',al,'N',N,'maxFE',maxFE);
            if i == 1 %SO通过全局变量result_t记录
                best(i,j,k) = result_t(end);
            end
            obj_min(i,j,k) = min(obj);
        end
    end
end
