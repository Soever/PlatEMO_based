function [Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f,GYbest,Xfood] = updateXbest(Xm,Xf,fitness_m,fitness_f,Xbest_m,Xbest_f,fitnessBest_m,fitnessBest_f)
    [Ybest1,gbest1] = min(fitness_m);
    [Ybest2,gbest2] = min(fitness_f);
    if Ybest1<fitnessBest_m
        Xbest_m = Xm(gbest1,:);
        fitnessBest_m=Ybest1;
    end
    if Ybest2<fitnessBest_f
        Xbest_f = Xf(gbest2,:);
        fitnessBest_f=Ybest2;
    end
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        Xfood=Xbest_m;
    else
        GYbest=fitnessBest_f;
        Xfood=Xbest_f;
    end
end
