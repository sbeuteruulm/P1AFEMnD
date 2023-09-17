function etaR = computeEtaRFull(x,coordinates,elements,element2neighbour,f,g,...
    D,volumes,grad)

nE = size(elements,1);
nD = size(coordinates,2);
%*** Compute grad u_h for each element
G = zeros(nE,nD);
for k = 1:nD+1
  G = G + grad{k}.*x(elements(:,k),ones(1,nD));
end
%*** Compute du_h/dn for each element and hyperface
dudn = zeros(nE,nD+1);
facets = flipud(nchoosek(1:nD+1,nD));
for k=1:nD+1
  %*** Compute jump terms
  J = zeros(nE,nD);
  idx = element2neighbour(:,k) > 0;
  J(idx,:) = (G(idx,:) - G(element2neighbour(idx,k),:))*D.A;
  %*** Index for Neumann boundary
  idx = element2neighbour(:,k) == -2;
  J(idx,:) = G(idx,:)*D.A;
  dudn(:,k) = -dot(grad{k},J,2)./sqrt(sum(grad{k}.^2,2));
  %*** Contributions on Neumann boundary
  if sum(idx) > 0
    [bary,wg] = quadnd(nD-1,5);
    neumann = elements(idx,facets(k,:));
    C = reshape(coordinates(neumann',:),nD,[]);
    L = zeros(size(neumann,1),1);
    for j = 1:size(bary,1)
      L = L + (g(reshape((bary(j,:)*C)',[],nD))-dudn(idx,k)).^2*wg(j);
    end
    dudn(idx,k) = sqrt(L);
  end
end
%*** Calculate volume terms
%*** Residual terms b dot grad u_h
bdu = G*D.b;
%*** Quadrature formula
[bary,wg] = quadnd(nD,2);
C = reshape(coordinates(elements',:),nD+1,[]);
fint = zeros(nE,1);
for j = 1:size(bary,1)
  fint = fint + (bdu + D.c*x(elements)*bary(j,:)' - f(reshape((bary(j,:)*C)',[],nD))).^2 * wg(j);
end
%*** Determine residual based error estimator
etaR = volumes.*((factorial(nD)*volumes).^(2/nD).*fint+nD*sum(dudn.^2,2));