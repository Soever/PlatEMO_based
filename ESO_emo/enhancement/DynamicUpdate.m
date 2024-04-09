function [C1,C2,C3] = DynamicUpdate(C1,C2,C3)
    T = length(C1) ;
    for l= 1:T
        C1(1,l)=C1(1,l)+0.1*sin((pi/2)*((rand).^4));%eq.(14)
        C2(1,l)=C2(1,l)+0.001*cos((pi/2)*((rand).^4));%eq.(15)
        C3(1,l)=C3(1,l)-2*sin(0.5*pi*(l/T)^4);%eq.(16)
    end
end

