function [x,energy,vol,G] = solvePDEFull(coordinates,elements, ...
                                     dirichlet,neumann,f,g,uD,D)
% Calculation of an approximated solution of the a general linear elliptic
% partial differential equation
%
% [x,energy,vol,G] = solvePDE(coordinates,elements, ...
% dirichlet,neumann,f,g,uD,D) calculates approximated solution x of PDE on
% a given mesh. The mesh is described by its coordinates given vertex wise
% in the rows of the array coordinates, the element array containing the
% vertex numbers for each elements in its rows and the boundary
% information, where dirichlet and neumann have in each row a hyperface of
% the boundary. The problem defining functions f, g, and uD are given as
% function handles. And the problem parameters are saved in object D which
% consists of the fields D.A, D.b, and D.c.
%
%Comments:
%   The function expects a non-empty array for dirichlet boundary
%   conditions. Neumann boundary conditions are optional for the solution
%   and can therefore be given as an empty array.
%
%Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
%Authors:
%   S. Beuter, S. Funken 18-10-22
nC = size(coordinates,1);
nD = size(coordinates,2);
%*** Compute gradients and relative volume of simplices |T|/|T_ref|
[vol,G] = volngrad(coordinates,elements);
%*** Assemble,   -div( A * Du ) + b * Du + c * u = f
S = sparse(nC,nC);
cc = D.c/((nD+1)*(nD+2)) * vol;
for j = 1 : nD+1
  volAxGj = vol.*( G{j} * D.A );  
  bbpcc = vol.*( G{j} * D.b/(nD+1) ) + cc;
  S = S + sparse(elements(:,j),elements(:,j),dot(G{j},volAxGj,2) + bbpcc + cc,nC,nC); 
  for k = j+1 : nD+1
    volGkGj = dot(G{k},volAxGj,2);
    S = S + sparse(elements(:,j),elements(:,k),volGkGj,nC,nC) ...
        + sparse(elements(:,k),elements(:,j),volGkGj,nC,nC);
    S = S + sparse(elements(:,k),elements(:,j),bbpcc,nC,nC); 
  end
  for k = 1:j-1
    S = S+sparse(elements(:,k),elements(:,j),bbpcc,nC,nC);
  end
end
%*** Prescribe values at Dirichlet nodes
x = zeros(nC,1);
dirichlet = unique(dirichlet);
x(dirichlet) = uD(coordinates(dirichlet,:));
%*** Assembly of volume force
[bary,wg] = quadnd(nD,2);
C = reshape(coordinates(elements',:),nD+1,[]);
L = zeros(size(elements,1),nD+1);
for k = 1:size(bary,1)
  L = L + (f(reshape((bary(k,:)*C)',[],nD)).*(vol*wg(k)))*bary(k,:);
end
b = accumarray(elements(:),L(:),[nC 1]) - S * x;
%*** Assembly of Neumann load
if size(neumann,1) 
  volNeu = volngrad(coordinates,neumann);
  [bary,wg] = quadnd(nD-1,1);
  C = reshape(coordinates(neumann',:),nD,[]);
  N = zeros(size(neumann,1),nD);
  for k = 1:size(wg,1)
    N = N + (g(reshape((bary(k,:)*C)',[],nD)).*(volNeu*wg(k)))*bary(k,:);
  end
  b = b + accumarray(neumann(:),N(:),[nC 1]);
end
%*** Computation of P1-FEM approximation
freenodes = setdiff(1:nC, dirichlet);
afun = @(x) S(freenodes,freenodes) * x;
D = diag(S(freenodes,freenodes));
prae = @(x) x./D;
x(freenodes) = pcg(@(x)afun(x),b(freenodes),1e-6*length(freenodes)^(-1/nD), ...
                  length(freenodes),@(x)prae(x) );
%*** Compute energy || grad(uh) ||^2 of discrete solution
energy = x'*S*x;



