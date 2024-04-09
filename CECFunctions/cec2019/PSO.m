%% Particle Swarm Optimization - PSO
function [gBestScore,gBest,cg_curve]=PSO(noP,Max_iteration,lb,ub,dim,fobj)

% Define the PSO's paramters
ub = ub.*ones(1,dim);
lb = lb.*ones(1,dim);

wMax=0.9;
wMin=0.2;
c1=2;
c2=2;
vMax = (ub - lb) * 0.2;
% vMax = 6.*ones(1,dim);
vMin  = -vMax;

% Initializations
iter=Max_iteration;
vel=zeros(noP,dim);
pBestScore=zeros(noP);
pBest=zeros(noP,dim);
gBest=zeros(1,dim);
cg_curve=zeros(1,iter);
vel=zeros(noP,dim);
pos=zeros(noP,dim);

%Initialization
for i=1:size(pos,1)
    for j=1:size(pos,2)
        pos(i,j)=(ub(j)-lb(j))*rand()+lb(j);
        vel(i,j)=rand();
    end
end

for i=1:noP
    pBestScore(i)=inf;
end

% Initialize gBestScore for a minimization problem
gBestScore=inf;


for t=1:iter
    
    for i=1:size(pos,1)
        % Return back the particles that go beyond the boundaries of the search
        % space
        Flag4ub=pos(i,:)>ub;
        Flag4lb=pos(i,:)<lb;
        pos(i,:)=(pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        %Calculate objective function for each particle
        fitness=fobj(pos(i,:));
        
        if(pBestScore(i)>fitness)
            pBestScore(i)=fitness;
            pBest(i,:)=pos(i,:);
        end
        if(gBestScore>fitness)
            gBestScore=fitness;
            gBest=pos(i,:);
        end
    end
    
    %Update the W of PSO
    w=wMax-t*((wMax-wMin)/iter);
    %Update the Velocity and Position of particles
    for i=1:size(pos,1)
        for j=1:size(pos,2)
            vel(i,j)=w*vel(i,j)+c1*rand()*(pBest(i,j)-pos(i,j))+c2*rand()*(gBest(j)-pos(i,j));
            
            if vel(i,j)>vMax(j)
                vel(i,j)=vMax(j);
            end
            if vel(i,j)<-vMax(j)
                vel(i,j)=-vMax(j);
            end
            pos(i,j)=pos(i,j)+vel(i,j);
        end
    end
    cg_curve(t)=gBestScore;
end


end