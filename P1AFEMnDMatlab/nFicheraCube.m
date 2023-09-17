function [coordinates,elements,boundary] = nFicheraCube(n,s,varargin)
[coordinates,elements] = nSubCube(n,2*s+1,varargin{:});
coordinates = -1+2*coordinates;
%*** Delete simplices with midpoints > 0 (i.e. each component is positive)
midpts = zeros(size(elements,1),n);
for j=1:n+1
  midpts = midpts + coordinates(elements(:,j),:);
end
midpts = midpts./(n+1);
idx = logical(sum(midpts<0,2));
elements = elements(idx,:);
nodes = unique(elements);
i2fi = zeros(size(coordinates,1),1);
i2fi(nodes) = 1:length(nodes);
coordinates = coordinates(nodes,:);
elements = i2fi(elements);
%*** Extract boundary
boundary=[];
for k = 1 : n+1
  bdry = elements(:,[1:k-1,k+1:n+1]);
  midpts = zeros(size(bdry,1),n);
  for j=1:n
    midpts = midpts + coordinates(bdry(:,j),:);
  end
  midpts = midpts./n;
  idx = (max(abs(midpts),[],2) >=1-1e-4) | (min(midpts,[],2)>=-1e-4);
  boundary = [boundary;bdry(idx,:)];
end


