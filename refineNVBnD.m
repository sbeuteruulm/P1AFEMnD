function [coordinates,elements,element2neighbour,level] = ...
                           refineNVBnD(coordinates,elements,element2neighbour,...
                           level,varargin)
% This function gives a framework for the refinement process.
%
% [coordinates,elements,element2neighbour,level] = ...
% refineNVBnD(coordinates,elements,element2neighbour,level,varargin)
% updates a mesh by the newest vertex bisection. The mesh is described by
% the arrays coordinates containing the coordinates for each vertex row
% wise, the array elements that describes the simplices by their vertex
% numbers, element2neighbour assembles for each elements all its neighbours
% in one row and vector level that gives the level of each simplex. These
% arrays get enlarged and updated before, they are returned. The cell array
% varargin assembles boundary arrays and as last element a vector with
% element numbers that are marked for refinement.
%
% Comments:
%   This function is called to perform the REFINE step of the adaptive FEM.
%   All marked elements are bisected by the newest vertex bisection. To
%   ensure that every marked element and neighbours sharing the refinement
%   edge are bisected exactly once, the refine function is called, if the
%   global level has not been updated before. The refine function
%   identifies all neighbours sharing a refinement edge and initiates the
%   Newest Vertex Bisection.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S. Beuter, S. Funken 18-10-22

global C E L E2N i2fi nC nE  
nB = length(varargin)-1;
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
    refine(e);  
  end
end
coordinates = C(1:nC,:); elements = E(1:nE,:);  level = L(1:nE); E2N = E2N(1:nE,:);
element2neighbour = E2N;
level = L;