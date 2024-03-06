function GCPD = cal_GCPD(m,popObj,lb,d,g)
    GCPD = 0 ;
    for k = 1:m
        GCPD = GCPD + ((popObj(k)-(lb(k)+g(k)*d(k)))/d(k))^2 ;
    end
    GCPD = sqrt(GCPD);
end

