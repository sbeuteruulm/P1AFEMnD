function varargout = extractBoundary(elements,element2neighbour,nB)
%*** extract boundary data from updated element2neighbour
nE = size(elements,1);
nN = size(elements,2);
facets = flipud(nchoosek(1:nN,nN-1)-1);
varargout = cell(nB,1);
for k = 1:nB
  [i,j,~] = find(element2neighbour == -k);
  idx = i * ones(1,nN-1) + facets(j,:) * nE;
  varargout{k} = reshape( elements(idx) , length(i) , []);
end
