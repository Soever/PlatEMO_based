function [outputArg1,outputArg2] = multi_resolution_grid(PopObj,N,div1,div2)

    %非支配解的个数
    NoP = size(PopObj,1);
    
    %% Calculate the grid location of each solution
    fmax = max(PopObj,[],1);
    fmin = min(PopObj,[],1);
    grid_lb = fmin - (fmax-fmin)/(div1*2) ;
    grid_ub = fmax + (fmax-fmin)/(div1*2) ;
    d    = (fmax-fmin)/div;
    fmin = repmat(fmin,NoP,1);
    d    = repmat(d,NoP,1);
    GLoc = floor((PopObj-fmin)./d);
    GLoc(GLoc>=div) = div - 1;
    GLoc(isnan(GLoc)) = 0;
    
    %% Detect the grid of each solution belongs to
    [~,~,Site] = unique(GLoc,'rows');

    %% Calculate the crowd degree of each grid
    CrowdG = hist(Site,1:max(Site));
    
    %% Roulette-wheel selection
    TheGrid = RouletteWheelSelection(N,CrowdG);
    REP     = zeros(1,N);
    for i = 1 : length(REP)
        InGrid = find(Site==TheGrid(i));
        Temp   = randi(length(InGrid));
        REP(i) = InGrid(Temp);
    end
end





