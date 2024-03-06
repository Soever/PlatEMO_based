function REP = multi_resolution_grid(PopObj,N,div1,div2)

    %非支配解的个数
    NoP = size(PopObj,1);
    m = size(PopObj,2);
    
    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);
    fmin = min(PopObj,[],1);
    grid_lb = fmin - (fmax-fmin)/(div1*2) ;
    grid_ub = fmax + (fmax-fmin)/(div1*2) ;
    d_coarse    = (grid_ub-grid_lb)/div1;   %得到每个粗分辨率格子的间距
    grid_lb_coarse = repmat(grid_lb,NoP,1); %为每个非支配解分配下限
    d_coarse    = repmat(d_coarse,NoP,1);          %为每个非支配解分配间距
    GLoc_coarse = floor((PopObj-grid_lb_coarse)./d_coarse); %计算每一个维度属于哪一个格子
    GLoc_coarse(GLoc_coarse>=div1) = div1 - 1;   
    GLoc_coarse(isnan(GLoc_coarse)) = 0;
    
    %% Detect the grid of each solution belongs to
    [~,~,Site] = unique(GLoc_coarse,'rows');%返回[唯一行的升序排列C，C在GLoc_coarse第一次出现的行号，GLoc_coarse行在C的对应行索引]

    %% 计算每个粗分辨率格子的拥挤度（存在非支配解的个数）
    CrowdG = hist(Site,1:max(Site)); %crowdg的下标是第i种格子，值是格子中解的个数
    max_CrowdG = max(CrowdG) ;
    max_CrowdG_idx = find(CrowdG==max_CrowdG) ;
   
    %% 找到第i种格子对应的个体
    if length(max_CrowdG_idx) == 1
        pop_index  = find(Site==max_CrowdG_idx);
        pop_num  = length(pop_index) ;
        pop_coarse = PopObj(pop_index,:);
        
        for j=1:length(pop_index)
            k = pop_index(j);
            GCPD(j)=cal_GCPD(m,PopObj(k,:),grid_lb_coarse(k,:),d_coarse(k,:),GLoc_coarse(k,:));
        end
        pop_delete = pop_index(find(GCPD==max(GCPD)));

    else
        for i= 1:length(max_CrowdG_idx)
            pop_index  = find(Site==max_CrowdG_idx(i));
            pop_num  = length(pop_index) ;
            pop_coarse = PopObj(pop_index,:);
            for j=1:length(pop_index)
                k = pop_index(j);
                GCPD(j)=cal_GCPD(m,PopObj(k,:),grid_lb_coarse(k,:),d_coarse(k,:),GLoc_coarse(k,:));
            end
            pop_delete = pop_index(find(GCPD==max(GCPD)));
        end
    end
    
    %% Roulette-wheel selection
    TheGrid = RouletteWheelSelection(N,CrowdG);
    REP     = zeros(1,N);
    for i = 1 : length(REP)
        InGrid = find(Site==TheGrid(i));
        Temp   = randi(length(InGrid));
        REP(i) = InGrid(Temp);
    end
end





