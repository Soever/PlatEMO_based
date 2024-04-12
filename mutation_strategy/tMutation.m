function newX_dec= tMutation(X_dec,df)
     t_randoms = trnd(df, size(X_dec, 1));
     newX_dec = Xm .* (1 + t_randoms);
end

