function [coordinates,elements] = nSubCube(n,s,varargin)
% Creates a triangulation of n-dimensional cube with edge length 1 which is
% split into n! x s^n Kuhn simplices in subcubes.
%
% [coordinates,elements] = nSubCube(n,s,varargin) creates a
% triangulation of a cube of edge length 1 split into s x s x s ... x s
% subcubes. Each subcube is split into Kuhn simplices. The subcubes are
% reflected, if varargin is 'tucker' and translated if varargin is not
% defined or 'freudenthal'. The function returns an array coordinates
% containing the vertex coordinates, and an array elements containing the
% simplices defined by their node numbers in each row.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S. Beuter, S. Funken 18-10-22

%*** Create coordinates
coordinates = fliplr(convertdec(n,s+2)/(s+1));
%*** Create elements
elem = KuhnTriangulation(n);
if s == 0 
  elements = elem;
  return
end
nE = size(elem,1);
A = fliplr(convertdec(n,s+1));
elements = zeros(((s+1)^n)*nE,n+1);
i2fi = [1,2];
for k=1:n-1
  i2fi = [i2fi,i2fi+(s+2)^k];
end

if nargin < 3 || strcmp( lower(varargin{1}) , 'freudenthal')
  pow = (s+2).^(0:n-1);
  for k=1:(s+1)^n
    tmp = pow*A(k,:)'+i2fi(elem);
    elements((k-1)*nE+(1:nE),:) = tmp;
  end
elseif strcmp(lower(varargin{1}),'tucker')
  map = ones(2^n+1,1)*i2fi;
  for k = 1:n
    jdx = [];
    for j=1:2^(n-k) 
      jdx = [jdx,(j-1)*2^k+(1:2^(k-1))];
    end
    idx = jdx + 2^(k-1);
    % exchange of colums idx and jdx, for each row idx
    tmp = map(idx,jdx);
    map(idx,jdx) = map(idx,idx);
    map(idx,idx) = tmp;
  end
  % 
  pow2 = 2.^(0:n-1);
  pow = (s+2).^(0:n-1);
  for k=1:(s+1)^n
    elements((k-1)*nE+(1:nE),:) = pow*A(k,:)' + reshape(...
                               map(1+pow2*mod(A(k,:)',2),elem),[],n+1);
    
  end
end
