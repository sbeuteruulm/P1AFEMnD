function [x,energy,vol,G] = solvePDEsymII(coordinates,elements, ...
                                     dirichlet,neumann,f,g,uD,D)
nC = size(coordinates,1);
nD = size(coordinates,2);
%*** Compute gradients and relative volume of simplices |T|/|T_ref|
[vol,G] = volngrad(coordinates,elements);
%*** Assemble,   -div( A * Du ) + b * Du + c * u = f
idx = repmat(1:nD+1,nD+1,1);
I = reshape(elements(:,idx(:)),[],1);idx = idx';
J = reshape(elements(:,idx(:)),[],1);
S = zeros(size(elements,1),(nD+1)^2);
% cc = D.c/((nD+1)*(nD+2)) * vol;
for j = 1 : nD+1
  volxGj = vol.*G{j};  
%   bbpcc = vol.*( G{j} * D.b/(nD+1) ) + cc;
  S(:,(j-1)*(nD+1)+j) = dot(G{j},volxGj,2); 
%   S = S+sparse(elements(:,j),elements(:,j),dot(G{j},volxGj,2),nC,nC)./2;
  for k = j+1 : nD+1
    S(:,(j-1)*(nD+1)+k) = dot(G{k},volxGj,2);
%     S = S+sparse(elements(:,j),elements(:,k),dot(G{k},volxGj,2),nC,nC);
    S(:,(k-1)*(nD+1)+j) = S(:,(j-1)*(nD+1)+k); % + bbpcc; 
    %S2 = S2+sparse(elements(:,k),elements(:,j),dot(G{k},volAxGj,2),nC,nC);
  end
%   for k = 1:j-1
%     S(:,(k-1)*(nD+1)+j) = S(:,(k-1)*(nD+1)+j);% + bbpcc; 
%   end
end
% S = S + S';
S = sparse(I,J,S(:));
% clear I J;
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



