% 假设你有一个 n*m 的数据矩阵 X
X = rand(10, 10);
% 步骤1: 创建一个KNN模型
% 示例数据，假设有一个 5*3 的矩阵，其中有 5 个个体，每个个体有 3 个特征
data =X; 

% 设置 k 值
k = 2;

% 初始化一个矩阵来存储每个个体最近的 k 个个体的索引
nearest_indices_all = zeros(size(data, 1), k+1);

% 遍历每个个体
for i = 1:size(data, 1)
    % 使用 knnsearch 函数找到每个个体最近的 k 个个体的索引
    nearest_indices = knnsearch(data, data(i, :), 'K', k+1); % +1 是因为包含了自身
    
    % % 排除自身
    % nearest_indices = nearest_indices(nearest_indices ~= i);
    % 
    % 将最近的 k 个个体的索引存储到矩阵中
    nearest_indices_all(i, :) = nearest_indices;
end

c = X(nearest_indices_all(1,:),:);
% 打印每个个体最近的 k 个个体的索引
disp(['每个个体最近的 ', num2str(k), ' 个个体的索引：']);
disp(nearest_indices_all);

% 现在 nearest_neighbors{i} 包含了第 i 个个体的近邻
