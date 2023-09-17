function boundary = nCubeBoundary(coordinates,elements)
nD = size(coordinates,2);
boundary=[];
for k = 1 : nD+1
  bdry = elements(:,[1:k-1,k+1:nD+1]);
  midpts = zeros(size(bdry,1),nD);
  for j=1:nD
    midpts = midpts + coordinates(bdry(:,j),:);
  end
  midpts = midpts./nD;
  boundary = [boundary;bdry( find(max(abs(midpts),[],2) >= 1-1e-4),:)];
  boundary = [boundary;bdry( find(min(abs(midpts),[],2) <=   1e-4),:)];
end