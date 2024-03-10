classdef MR_PSO < ALGORITHM
   methods
        function main(Algorithm,Problem)
            %% Parameter setting
            div = Algorithm.ParameterSet(10);

            %% Generate random population
            Population = Problem.Initialization();
            Archive    = UpdateArchive(Population,Problem.N,div);
            Pbest      = Population;

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                REP        = multi_resolution_grid(Archive.objs,Problem.N,div);
                Population = OperatorPSO(Problem,Population,Pbest,Archive(REP));
                Archive    = UpdateArchive([Archive,Population],Problem.N,div);
                Pbest      = UpdatePbest(Pbest,Population);
            end
        end
    end
end
