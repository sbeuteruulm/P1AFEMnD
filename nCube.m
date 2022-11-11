function [coordinates,elements,boundary] = nCube(n,s,varargin)
% Creats a triangulation with explicit boundary numbers of n! x s^n Kuhn 
% simplices in a cube
%
% [coordinates,elements,boundary] = nCube(n,s,varargin) creats a
% triangulation of a cube of edge length 1 split into s x s x s ... x s
% subcubes. Each subcube is split into Kuhn simplices. The subcubes are
% reflected, if varargin is 'tucker' and translated if varargin is not
% defined or 'freudenthal'. The function returns an array coordinates
% containing the vertex coordinates, an array elements containing the
% simplices defined by their node numbers in each row and an array boundary
% giving the boundary nodes.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S.Beuter, S. Funken 18-10-22

[coordinates,elements] = nSubCube(n,s,varargin{:});
boundary=[];
for k = 1 : n+1
  bdry = elements(:,[1:k-1,k+1:n+1]);
  midpts = zeros(size(bdry,1),n);
  for j=1:n
    midpts = midpts + coordinates(bdry(:,j),:);
  end
  midpts = midpts./n;
  idx = (max(abs(midpts),[],2) >=1-1e-4) | (min(midpts,[],2)<=1e-4);
  boundary = [boundary;bdry(idx,:)];
end


