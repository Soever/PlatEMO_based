population = rand(50, 10);
% 设置 K-means++ 算法的选项
opts = statset('Display','final');

% 使用 K-means++ 算法对种群进行聚类
k =2 ;
[idx, centroids] = kmeans(population, k, 'Options',opts, 'Start','plus');

% idx 表示每个个体所属的类别编号，centroids 表示每个类别的中心点坐标

% 可视化聚类结果
% 例如，对于二维决策变量，可以使用散点图表示不同类别的个体
scatter(population(:,1), population(:,2), 10, idx, 'filled');
hold on;
scatter(centroids(:,1), centroids(:,2), 50, (1:k)', 'filled', 'MarkerEdgeColor', 'k');
legend('个体', '中心点');
title('K-means++ 聚类结果');
xlabel('决策变量1');
ylabel('决策变量2');