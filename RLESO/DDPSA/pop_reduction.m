function newX = pop_reduction(X,fitness,lamda)
    [~,gbest1] = min(fitness);
    X_best = X(gbest,:);
    %% 删除最佳X
    indices = true(size(fitness));indices(gbest1) = false;
    fitness = fitness(indices);
    X = X(indices,:);
    %% 基于fitness计算生存度
    maxFit  = max(fitness) ;minFit = min(fitness) ;
    val_fit = (maxFit-fitness) ./ (maxFit-minFit+eps) ; % 适应度越小 生存度越大 个体越好
    %% 基于距离计算生存度
    dis =  sqrt(sum((X-X_best).^2,2));
    minDis = min(dis) ;maxDis = max(dis) ;
    val_div = (dis-minDis) ./  (maxDis-minDis+eps);% 距离越大 生存度越大 个体越好
    %% 总生存度
    val =lamda*val_fit+(1-lamda)*val_div;
    
    [sortedValues, sortedIndices] = maxk(val , 3);
    best3X = X(sortedIndices,:);
    best3X_fitness = fitness(sortedIndices);
    w = (sum(best3X_fitness)-best3X_fitness)/2*sum(best3X_fitness) ;
    C =sum(w.*best3X,1) ;
    r1 = rand ;
    newX = (1-r1)*X_best+r1*C;
end

