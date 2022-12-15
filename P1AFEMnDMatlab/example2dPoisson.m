clear
close all
clc

% Define 2d mesh with L-shaped domain ([-1,1]x[-1,1])/([0,1]x[0,1])
coordinates = [-1 -1;
                0 -1;
               -1  0;
                0  0;
                1  0;
               -1  1;
                0  1;
                1  1];
            
elements = [1 2 4;
            1 3 4;
            6 3 4;
            6 7 4;
            4 7 8;
            4 5 8];
        
dirichlet = [2 1;
             1 3;
             3 6;
             6 7;
             7 8;
             8 5];
         
neumann = [2 4;
           4 5];
       
level = zeros(size(elements,1),1);

% Define problem specification
D.A  = eye(2);
D.b  = zeros(2,1);    
D.c  = 0; 
            
% Define right hand side
f  = @(x)  ones(size(x,1),1);
g  = @(x)  zeros(size(x,1),1);
uD = @(x)  zeros(size(x,1),1);

% Define parameters for stopping criterion
nEmax = 10000;
eta = 0.7;

% Run FEM
[x,energy,coordinates,elements,bdry] = adaptiveAlgorithm( ...
                   coordinates,elements,level,{dirichlet,neumann},f,g,uD,D,nEmax,eta);
      
% Plot solution
trisurf(elements,coordinates(:,1),coordinates(:,2),x,'facecolor','interp')