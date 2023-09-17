function [coordinates,elements,boundary] = nCube(n,s,varargin)
% Creates a triangulation with explicit boundary numbers of n! x s^n Kuhn 
% simplices in a cube
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


