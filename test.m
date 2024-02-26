maxIteration = 500 ;
N = 30 ;
D = 30 ;
maxFE = maxIteration*N ;
platemo('problem',@CEC2017_F4,'algorithm',@SnakeOptimizer,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2017_F4,'algorithm',@GA,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2017_F4,'algorithm',@PSO,'N',N,'D',D,'maxFE',15000);
platemo('problem',@CEC2017_F4,'algorithm',@SHADE,'N',N,'D',D,'maxFE',15000);