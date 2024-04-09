function Temp = ConvexLensImaging(k,Xbest,ub,ub1,lb,lb1)
    dim = length(Xbest);
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

