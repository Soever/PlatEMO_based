clc
clear

e = exp(1);

[Dec,obj,Con]=platemo('objFcn',@ObjFcn,'algorithm',@GA,'N',100,'D',5,'maxFE',10000,'lower',[0,0,1/e,0,0],'upper',[1,1,1,10,10]);


%% 定义评价函数
function f = ObjFcn(x)
    al = {@SnakeOptimizer,x(1),x(2),x(3),x(4),x(5)} ;
    res = zeros(5) ;
    parfor i = 1:5
        [~,obj1,~]=platemo('problem',@CEC2008_F1,'algorithm',al,'N',100,'maxFE',10000);
        [~,obj4,~]=platemo('problem',@CEC2008_F4,'algorithm',al,'N',100,'maxFE',10000);
        [~,obj7,~]=platemo('problem',@CEC2008_F7,'algorithm',al,'N',100,'maxFE',10000);
        res(i) = (min(obj1)+min(obj4)+min(obj7))/3;
    end
    f = mean(res) ;
end
