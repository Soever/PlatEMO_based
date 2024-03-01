meanMatrix  = zeros(10,4);
varMatrix = zeros(10,4);

for ii = 1:4
    for jj = 1:10
        slice = squeeze(obj_min(ii, jj, :)); % 获取第ii个i和第jj个j对应的k维数据
        meanMatrix(jj, ii) = mean(slice); % 计算平均值
        varMatrix(jj, ii) = var(slice); % 计算方差
    end
end