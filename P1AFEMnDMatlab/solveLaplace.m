function [u,energy] = solveLaplace(coordinates,elements,dirichlet,neumann,f,g,uD)
% Initialisation
nD = size(coordinates,2); 
% Assembly
A = sparse(size(coordinates,1),size(coordinates,1));
for j = 1:size(elements,1)
  A(elements(j,:),elements(j,:)) = A(elements(j,:),elements(j,:)) ...
      + stima(coordinates(elements(j,:),:));
end
% Volume Forces
b = zeros(size(coordinates,1),1);
for j = 1:size(elements,1)
  b(elements(j,:)) = b(elements(j,:)) + ...
      abs(det([ones(1,nD+1);coordinates(elements(j,:),:)'])) ... 
      * f(sum(coordinates(elements(j,:),:))/(nD+1)) / factorial(nD+1);
end
% Neumann conditions
for j = 1 : size(neumann,1)
  tmp = coordinates(neumann(j,2:end),:) ...
      - coordinates(neumann(j,ones(1,nD-1)),:);
  b(neumann(j,:)) = b(neumann(j,:)) + sqrt(det(tmp*tmp')) ...
      * g(sum(coordinates(neumann(j,:),:))/nD)/factorial(nD);
end
% Dirichlet conditions 
dirichlet = unique(dirichlet);
FreeNodes = setdiff(1:size(coordinates,1),dirichlet);
u = zeros(size(coordinates,1),1);
u(dirichlet) = uD(coordinates(dirichlet,:));
b = b - A * u;
% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);
energy = u' * A * u;
