clear
close all
clc

% Define settings
    % problem size ('small', 'medium', 'large')
    Psize = 'medium';
    % adaptivity ('low', 'high', 'uniform')
    adaptivity = 'high';
    % problem dimesnision (a natural number)
    d = 3;
    % problem type ('linear', 'quadratic' (in the first component))
    pType = 'quadratic';
    
% Define stopping criterion and adaptivity parameter according to problem
% setting
switch true
    case strcmp(Psize,'small') && strcmp(adaptivity,'low')
        nEmax = 1000;
        eta = 0.7;
    case strcmp(Psize,'small') && strcmp(adaptivity,'high')
        nEmax = 1000;
        eta = 0.2;
    case strcmp(Psize,'small') && strcmp(adaptivity,'uniform')
        nEmax = 1000;
        eta = 1;
    case strcmp(Psize,'medium') && strcmp(adaptivity,'low')
        nEmax = 50000;
        eta = 0.7;
    case strcmp(Psize,'medium') && strcmp(adaptivity,'high')
        nEmax = 50000;
        eta = 0.2;
    case strcmp(Psize,'medium') && strcmp(adaptivity,'uniform')
        nEmax = 50000;
        eta = 1;
    case strcmp(Psize,'large') && strcmp(adaptivity,'low')
        nEmax = 1000000;
        eta = 0.7;
    case strcmp(Psize,'large') && strcmp(adaptivity,'high')
        nEmax = 1000000;
        eta = 0.2;
    case strcmp(Psize,'large') && strcmp(adaptivity,'uniform')
        nEmax = 1000000;
        eta = 1;
    otherwise
        disp('No possible setting choice.')
        return
end

% Define mesh dependent on dimension d and subcube size 1/(s+1). The domain
% is a Fichera cube [-1,1]^d/[0,1]^d. Here the whole boundary is formed of
% dirichlet boundary, while the neumann boundary is empty.
s=1;
[coordinates,elements,dirichlet] = nFicheraCube(d,s,'tucker');
neumann = [];

level = zeros(size(elements,1),1);

% Define problem specification
D.A  = ones(d)+d*eye(d);
D.b  = zeros(d,1);
D.c  = 0;

% Define right hand side
if strcmp(pType,'linear')
    f  = @(x)  zeros(size(x,1),1);
    g  = @(x)  zeros(size(x,1),1);
    uD = @(x)  x(:,1);
elseif strcmp(pType,'quadratic')
    f  = @(x)  -(2*d+2)*ones(size(x,1),1);
    g  = @(x)  zeros(size(x,1),1);
    uD = @(x)  x(:,1).^2 + x(:,2);
else
    disp('problem type not supported.')
end

% Run FEM
[x,energy,coordinates,elements,bdry] = adaptiveAlgorithm( ...
                   coordinates,elements,level,{dirichlet,neumann},f,g,uD,D,nEmax,eta);

% Relative error
if strcmp(pType,'linear')
    [maxi,idx] = max(abs(x-coordinates(:,1)));
    if coordinates(idx,1) ~= 0
        disp(['Maximum relative error ',num2str(maxi/abs(coordinates(idx,1)))])
    end
else
    [maxi,idx] = max(abs(x-(coordinates(:,1).^2+coordinates(:,2))));
    if (coordinates(idx,1).^2+coordinates(idx,2)) ~= 0
        disp(['Maximum relative error ',num2str(maxi/abs(coordinates(idx,1).^2+coordinates(idx,2)))])
    end
end

% Plot solution
if size(coordinates,1) < 1e5
    if d == 2
        trisurf(elements,coordinates(:,1),coordinates(:,2),x,'facecolor','interp')
    elseif d == 3
        trisurf(bdry{1},coordinates(:,1),coordinates(:,2), ...
                coordinates(:,3),x,'facecolor','interp')
    end
end