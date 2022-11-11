function [x,energy,coordinates,elements,bdry] = adaptiveAlgorithm( ...
                   coordinates,elements,level,bdry,f,g,uD,D,nEmax,rho)
% Framework function for the realization of the adaptive finite element
% method in arbitrary dimensions
%
% [x,energy,coordinates,elements,bdry] = adaptiveAlgorithm( ...
%                   coordinates,elements,level,bdry,f,g,uD,D,nEmax,rho)
% Realization of the loop of the adaptive FEM. The mesh data is given as
% input arrays, where coordinates contains the vertex coordinates, elements
% contains the simplices given by their vertex numbers in each row, level
% is a vector with the refinement level of each element and bdry is a cell
% array giving the nodes of each boundary type. The model parameters are
% given by the functionhandles f,g,uD and the object D for the model
% parameters. The maximal number of elements nEmax and the refinement
% parameter rho define the stopping criterion and the adaptivity of the
% mesh. As output the mesh is updated and is return in the same format as
% the input mesh. Additionally, the solution x is given as a vector.
%
% Comments:
%   The function wraps the four steps SOLVE, ESTIMATE, MARK and REFINE of
%   the adaptive FEM in a loop. For the steps SOLVE, ESTIMATE and REFINE
%   functions are called to realize the respective steps. The algorithm
%   stops when a given number of elements is passed.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S. Beuter, S. Funken 18-10-22

nB = length(bdry);
%*** Initiate element2neighbour
[~,element2hyperfaces,bdry2hyperfaces{1:nB}] = provideHyperFaceData(elements,bdry{:});
element2neighbour = provideNeighbourData(element2hyperfaces,bdry2hyperfaces{:});         
while 1
  %*** Compute discrete solution
  [x,energy,vol,G] = solvePDE(coordinates,elements,bdry{:},f,g,uD,D);
  %*** Compute refinement indicators
  indicators = computeEtaR(x,coordinates,elements,element2neighbour,f,g,vol,G);
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
          refineNVBnD(coordinates,elements,level,bdry{:},marked);
  [bdry{:}] = extractBoundary(elements,element2neighbour,nB);
end
