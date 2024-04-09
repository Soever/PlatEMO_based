function  TempX= convex_lens_imaging(ub,lb,k,Xbest_dec,Problem)

    TempX_dec = (ub + lb)./2 + (ub + lb)./(2*k) - Xbest_dec./k;
    TempX_dec = Problem.CalDec(TempX_dec) ;
    TempX =Problem.Evaluation(TempX_dec) ;
    
end

