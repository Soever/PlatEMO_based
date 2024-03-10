function REP = multi_resolution_grid(PopObj,N,div1,div2)
    %非支配解的个数
    NoP = size(PopObj,1);
    m = size(PopObj,2);
    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);fmin = min(PopObj,[],1);
    lb = fmin - (fmax-fmin)/(div1*2) ;ub = fmax + (fmax-fmin)/(div1*2) ;
    d = (ub-lb)/div1;   %得到每个粗分辨率格子的间距
    lb= repmat(lb,NoP,1); %为每个非支配解分配下限
    d  = repmat(d,NoP,1);          %为每个非支配解分配间距
    %计算每一个维度属于哪一个格子，GLOC为Nop*m，value为坐标值
    GLoc_coarse = floor((PopObj-lb)./d); GLoc_coarse(GLoc_coarse>=div1) = div1 - 1; GLoc_coarse(isnan(GLoc_coarse)) = 0;
    %% Detect the grid of each solution belongs to
    [~,~,Site] = unique(GLoc_coarse,'rows');%返回[唯一行的升序排列C，C在GLoc_coarse第一次出现的行号，GLoc_coarse行在C的对应行索引]

    %% 计算每个粗分辨率格子的拥挤度（存在非支配解的个数）
    CrowdG = hist(Site,1:max(Site)); %crowdg的下标是第i种格子，值是格子中解的个数
    max_value = max(CrowdG) ;
    if max_value>1 
        Spies = find(CrowdG==max_value) ;%得到最大拥挤度的格子种类，value是site的第几种格子
        NoG =  length(Spies) ;
        while NoG ~= 1
            eveness = zeros(NoG,1) ;
            %pop_coarse = {} ;
            Site_new = zeros(size(Site));
            Site_fine = {} ;
            pop_index = {} ;
            for i = 1:NoG
                pop_index{i}  = find(Site==Spies(i));
                pop_coarse = PopObj(pop_index{i} ,:);
                [eveness(i),Site_fine{i}] = calculateEveness(div2,pop_coarse);
     
            end
            min_eveness = min(eveness) ;
            Spies = find(eveness==min_eveness) ;
            NoG =  length(Spies) ;
            if NoG ~= 1
                Spices_all = 0 ;
                for i = 1 :NoG
                    for j= 1:length(Site_fine{i})
                        Site_new(pop_index{i}(j)) = Spices_all+Site_fine{i}(j);
                    end 
                    Spices_all = Spices_all+max(Site_fine{i}) ;
                end
                Site = Site_new ;
                NoG = Spices_all;
                Spies = 1: Spices_all;
                CrowdG   = histcounts(Site, 0.5:1:Spices_all+0.5);
            end
            
        end
    
        pop_index  = find(Site==Spies);
        pop_num  = length(pop_index) ;
        pop_coarse = PopObj(pop_index,:);
        
        for j=1:length(pop_index)
            k = pop_index(j);
            GCPD(j)=cal_GCPD(m,PopObj(k,:),lb(k,:),d(k,:),GLoc_coarse(k,:));
        end
        pop_delete = pop_index(find(GCPD==max(GCPD)));
        index = [1:pop_index-1, pop_index+1:NoP];
        PopObj = PopObj(index,:);
    end
    REP = REPSelection(PopObj,N,div1);

    
    

end


function GCPD = cal_GCPD(m,popObj,lb,d,g)
    GCPD = 0 ;
    for k = 1:m
        GCPD = GCPD + ((popObj(k)-(lb(k)+g(k)*d(k)))/d(k))^2 ;
    end
    GCPD = sqrt(GCPD);
end

function [eveness,Site] = calculateEveness(div2,PopObj)
    pop_num = size(PopObj,1);
    m = size(PopObj,2);


    ub = max(PopObj,[],1);lb = min(PopObj,[],1);
    d = (ub-lb)/div2;
    lb = repmat(lb,pop_num,1); %为每个非支配解分配下限
    d  = repmat(d,pop_num,1);          %为每个非支配解分配间距
    GLoc = floor((PopObj - lb)./d); %计算每一个维度属于哪一个格子
    GLoc(GLoc>=div2) = div2 - 1;   
    GLoc(isnan(GLoc)) = 0;

    [~,~,Site] = unique(GLoc,'rows'); %计算精细分辨率格子们哪些有解
    Crowd = hist(Site,1:max(Site));%计算每个精细分辨率格子的拥挤度
     
    x = Crowd ./ pop_num ;
    eveness = - sum(x .* log(x))/log(div2^m) ;

end




