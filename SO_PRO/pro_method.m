function Xnew = pro_method(N,D,pop,x_BEST,Tau,Schedule)
    
    ResetZero = zeros(1,D);
    Xnew = pop.decs ;
    [~, sortedIndx] = sort([pop.objs]);
    Selection_rate = exp(-(1-Tau));
    for i=1:N
         
        k = N;
        if i < N
            k = sortedIndx(randi([i+1 N]));
        end
        [~,Candid_Behavs] = sort(Schedule(i,1:D),'descend');
        Landa = ceil(D*rand*Selection_rate);
        Selected_behaviors = Candid_Behavs(1:Landa);%
        if rand <  0.5 %(0.1 + 0.9 * (1-Tau))
          Stimulation = ResetZero;
          Stimulation(Selected_behaviors) = ( x_BEST.dec(Selected_behaviors) - pop(i).dec(Selected_behaviors));    
        else
          Stimulation = ResetZero;
          Stimulation(Selected_behaviors) = ( pop(i).dec(Selected_behaviors) - pop(k).dec(Selected_behaviors));            
        end
        SF = Tau + rand * (mean(Schedule(i,Selected_behaviors) )/max(abs(Schedule(i,1:D))));
        Xnew(i,Selected_behaviors) = pop(i).dec(Selected_behaviors) + SF .* Stimulation(Selected_behaviors); 
    end

end

