classdef LA_PSO < ALGORITHM
    properties
        arch = {}; % Archive to store solution pairs
        arch_size = 100; % Maximum size of the archive
        lp = 0.4; % Threshold for learning-aided evolution
    end
    methods
        function main(Algorithm,Problem)
            ANN = feedforwardnet(10); % 创建一个具有两个隐藏层，每层10个神经元的前馈网络
            ANN.trainParam.showWindow = false;
            % ANN = configure(ANN,Problem.D,Problem.D); % 配置网络输入和目标
            
            %% Parameter setting
            
            [W,alpha,beta] = Algorithm.ParameterSet(0.4,0.5,0.9);
            %% Generate random population
            Population = Problem.Initialization();
            Pbest      = Population;
            [~,best]   = min(FitnessSingle(Pbest));
            Gbest      = Pbest(best);
            %% Optimization
            t = 0 ;
            while Algorithm.NotTerminated(Population)
                oldPopulation = Population ;
                t = t+1 ;
                r = rand() ;
                if r<Algorithm.lp & t>1 & size(Algorithm.arch,1)>0
                     Population = Opertor_LearningEvolution(Problem, Population.decs, ANN,Pbest.decs,alpha,beta );
                else
                     Population = OperatorPSO(Problem,oldPopulation,Pbest,Gbest,W);
                end
                replace        = FitnessSingle(Pbest) > FitnessSingle(Population);
                Pbest(replace) = Population(replace);
                [~,best]       = min(FitnessSingle(Pbest));
                Gbest          = Pbest(best);
                % SEP收集及存档更新
               
                for i = 1 : numel(Population)
                    if oldPopulation(i).obj < Population(i).obj
                        Algorithm.arch{end+1,1} =scale(oldPopulation(i).decs,Problem.lower,Problem.upper) ;
                        Algorithm.arch{end,2} = scale(Population(i).decs,Problem.lower,Problem.upper);
                    
                    end
                end
                % 如果需要，修剪存档
                if size(Algorithm.arch,1) > Algorithm.arch_size
                    Algorithm.arch = Algorithm.arch(end-Algorithm.arch_size+1:end,:);
                end
                % 训练ANN
                ANN= trainANN(ANN,Algorithm.arch);
            end
        end
    end
end

function ANN = trainANN(ANN,arch)
    if size(arch,1) > 0
        inputs = cell2mat(arch(:, 1));
        outputs = cell2mat(arch(:, 2));
        
        % trainOpts = trainingOptions('Plots','none', 'Verbose',false);
    
        ANN= train(ANN,inputs',outputs');
    end

end


function Population = Opertor_LearningEvolution(Problem, decs, ANN,Pbest_decs,alpha,beta)
    [n,d] = size(Pbest_decs);
    %learning-aided mutation
    pairs = generateRandomPairs(n) ;
    scale_Pbest_decs= scale(Pbest_decs,Problem.lower,Problem.upper);
    ind_b = pairs(:,1);ind_c = pairs(:,2);
    Xa = ANN(scale_Pbest_decs');
    Xa = Xa';
    Xa = scale_reverse(Xa,Problem.lower,Problem.upper);
    Xb = decs(ind_b,:);
    Xc = decs(ind_c,:);          
    newX = Xa +alpha*(Xb-Xc);
    %learning-aidid crossover
    r = rand(n,d);
    newX(r<=beta) =  Pbest_decs(r<=beta);

    Population=Problem.Evaluation(newX) ;

end
function pairs = generateRandomPairs(n)
    pairs = zeros(n, 2); % 初始化结果矩阵
    for i = 1:n
        validChoices = setdiff(1:n, i); % 移除当前索引，保留有效的选择
        idx = randperm(length(validChoices), 2); % 随机选取两个不同的索引
        pairs(i, :) = validChoices(idx); % 填充选取的值
    end
end
function scaled_P = scale(P,lower,upper)
% 首先将 lower 和 upper 扩展成 n x m 维矩阵
    lower_mat = repmat(lower, size(P, 1), 1);
    upper_mat = repmat(upper, size(P, 1), 1);
    % 然后执行缩放
    scaled_P = (P - lower_mat) ./ (upper_mat - lower_mat);
end
function scaled_P = scale_reverse(P,lower,upper)
% 首先将 lower 和 upper 扩展成 n x m 维矩阵
    lower_mat = repmat(lower, size(P, 1), 1);
    upper_mat = repmat(upper, size(P, 1), 1);
    % 然后执行缩放
    scaled_P = P  .* (upper_mat - lower_mat)+lower_mat;
end
