%% cec2019

function [lb,ub,dim,fobj] = Get_Functions_cec2019(F)

switch F
    case 1
        dim=9;
        lb=-8192*ones(1,dim);
        ub=8192*ones(1,dim);
        fobj = @(x) cec19_func(x',1);
        
    case 2
        dim=16;
        lb=-16384*ones(1,dim);
        ub=16384*ones(1,dim);
        fobj = @(x) cec19_func(x',2);
        
    case 3
        dim=18;
        lb=-4*ones(1,dim);
        ub=4*ones(1,dim);
        fobj = @(x) cec19_func(x',3);
        
        
    case 4
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec19_func(x',4);
        
    case 5
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec19_func(x',5);
        
        
    case 6
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec19_func(x',6);
        
    case 7
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec19_func(x',7);
        dim=10;
    case 8
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec19_func(x',8);
        dim=10;
    case 9
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec19_func(x',9);
        
    case 10
        dim=10;
        lb=-100*ones(1,dim);
        ub=100*ones(1,dim);
        fobj = @(x) cec19_func(x',10);
        
end

end
