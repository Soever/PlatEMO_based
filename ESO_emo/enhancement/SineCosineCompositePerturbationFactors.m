function [t1,t2] = SineCosineCompositePerturbationFactors(t1,t2,a)
    T = length(t1);
    for l= 1:T
        t1(1,l)=t1(1,l)+(0.0001*(sin(a*4*pi*l)+cos(a*6*pi*l)))*exp((pi/100)*(0.25*(T-l)));
        t2(1,l)=t2(1,l)+(0.0001*(cos(a*4*pi*l)+sin(a*6*pi*l)))*exp((pi/100)*(0.25*(T-l)));
    end
end

