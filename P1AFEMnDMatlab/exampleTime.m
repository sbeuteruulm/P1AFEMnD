clear
close all
clc

% Define settings
eta = 1; % uniform refinement
dmax = 4; % maximum dimension
d = 2:dmax; % dimensions to be calculated
reps = 4; % repeats for time measurements (skip first measurement)
nbRefinements = 10; % number of refinement steps

% allocate arrays
timeTmp = cell(length(d));
nbElements = zeros(nbRefinements+1,length(d));
tic
for i = d
    
    % Define problem specification
    D.A  = ones(i)+eye(i);
    D.b  = ones(i,1);    
    D.c  = 1; 

    % Define right hand side
    f  = @(x)  ones(size(x,1),1);
    g  = @(x)  zeros(size(x,1),1);
    uD = @(x)  zeros(size(x,1),1);
    
    timeTmp{i} = zeros(nbRefinements+1,reps);
    for j = 0:nbRefinements
        for k = 1:reps
            % Define mesh dependent on dimension d and subcube size 1/(s+1).
            % The domain is a cube [0,1]^d. Here the whole boundary is
            % formed of dirichlet boundary, while the neumann boundary is
            % empty.
            s=1;
            [coordinates,elements,dirichlet] = nCube(i,s,'tucker');
            neumann = [];

            level = zeros(size(elements,1),1);

            % Run FEM
            timeStart = toc;
            [x,energy,coordinates,elements,bdry] = adaptiveAlgorithm( ...
                               coordinates,elements,level,{dirichlet,neumann},...
                               f,g,uD,D,size(elements,1)*2^j,eta);
            timeTmp{i-1}(j+1,k) = toc-timeStart;  
        end
        nbElements(j+1,i-1) = size(elements,1);
    end
end

time = zeros(nbRefinements+1,length(d));
for i = 1:length(d)
    time(:,i) = mean(timeTmp{i}(:,2:end),2);
end

loglog(nbElements,time,'*-')
title('Time mesurements for adaptive FEM')
xlabel('number of elements')
ylabel('time in seconds')
set(gca,'FontSize',20)
grid on