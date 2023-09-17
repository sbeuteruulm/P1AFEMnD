function [x,energy,coordinates,elements,bdry] = adaptiveAlgorithm( ...
                   coordinates,elements,level,bdry,f,g,uD,D,nEmax,rho)

nB = length(bdry);
%*** Initiate element2neighbour
[~,element2hyperfaces,bdry2hyperfaces{1:nB}] = provideHyperFaceData(elements,bdry{:});
element2neighbour = provideNeighbourData(element2hyperfaces,bdry2hyperfaces{:}); 
while 1
  %*** Compute discrete solution
  [x,energy,vol,G] = solvePDEFull(coordinates,elements,bdry{:},f,g,uD,D);
  %*** Compute refinement indicators
  indicators = computeEtaRFull(x,coordinates,elements,element2neighbour,f,g,D,vol,G);
  
  %*** Stopping criterion
  if size(elements,1) >= nEmax
    break
  end
  %*** Mark elements for refinement
  [indicators,idx] = sort(indicators,'descend');
  sumeta = cumsum(indicators);
  ell = find( sumeta >= sumeta(end) * rho,1);
  marked = idx(1:ell);
  %*** Refine mesh 
  [coordinates,elements,element2neighbour,level] = ...
          refineNVBnD(coordinates,elements,element2neighbour,level,bdry{:},marked);
  [bdry{:}] = extractBoundary(elements,element2neighbour,nB);
end

