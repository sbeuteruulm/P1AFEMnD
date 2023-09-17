function [C,E,E2N,L] = refineNVBnD(coordinates,elements,element2neighbour,...
                           level,varargin)                
marked = varargin{end};

E2N = element2neighbour;
nE = size(elements,1); nC = size(coordinates,1);
C = coordinates; E = elements; L = level; 
%*** allocate memory for refined mesh
nEmax = ceil(2.2 * nE);
E(nEmax,1) = 0;  E2N(nEmax,1) = 0; L(nEmax,1) = 0; C(2*nC,1) = 0;
%*** refine mesh by bisection
lev = L;
for e = reshape(marked,1,[])
  if lev(e) == L(e)
    [~,C,E,L,E2N,nC,nE] = refine(e,C,E,L,E2N,nC,nE);  
  end
end
C = C(1:nC,:); E = E(1:nE,:);  L = L(1:nE); E2N = E2N(1:nE,:);