function Temp = ConvexLensImaging(Xbest,ub,ub1,lb,lb1)
    k = 10*(1-2*(t/T)^2);% scaling factor eq.(19)
    Temp = (ub + lb)./2 + (ub + lb)./(2*k) - Xbest./k; %eq.(20)
    for c = 1:dim
        if(Temp(1,c)>ub)
            Temp(1,c) =ub1(c);
        end
        if(Temp(1,c)<lb)
            Temp(1,c)  =lb1(c);
        end
    end
end

