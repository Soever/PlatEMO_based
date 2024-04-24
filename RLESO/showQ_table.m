function dataCell = showQ_table(A)
[n,m,~] = size(A) ;

% 创建一个大的图形窗口，并且根据需要调整大小
figure('Units', 'normalized', 'Position', [0, 0, 1, 1]); % 全屏
grayColor = [0.8, 0.8, 0.8];
% 对于每个元素位置
for i = 1:n
    for j = 1:m
        % 计算当前子图的位置
        subplot(n, m, (i-1)*m + j);
        
        % 取出当前位置的4个数值
        data = reshape(A(i, j, :), [2, 2]);
        
        % 设置当前坐标区域的边界，使得我们有空间放置数字
        set(gca, 'XLim', [0.5, 2.5], 'YLim', [0.5, 2.5]);
        
        % 关闭坐标轴
        axis off;
        
        % 遍历每个值，并在相应位置显示数值
        for k = 1:2
            for l = 1:2
                if data(k, l) ~= 0
                    patch([l-0.5, l+0.5, l+0.5, l-0.5], [k-0.5, k-0.5, k+0.5, k+0.5], grayColor, 'EdgeColor', 'none');
                end
                text(l, k, num2str(data(k, l), '%0.2f'), ...
                     'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle');
            end
        end
    end
end


end

