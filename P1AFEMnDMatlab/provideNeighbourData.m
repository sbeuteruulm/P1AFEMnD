function [element2neighbour,hyperface2elements] = ...
                          provideNeighbourData(element2hyperfaces,varargin)
nE = size(element2hyperfaces,1);
nN = size(element2hyperfaces,2);

hyperface2elements = zeros(max(element2hyperfaces(:)),2);
ind = sortrows([element2hyperfaces(:),repmat((1:nE)',nN,1)]);
flag = [false;ind(1:end-1,1) == ind(2:end,1)];
hyperface2elements(ind(~flag,1),1) = ind(~flag,2);
hyperface2elements(ind(flag,1),2) = ind(flag,2);

for j = 1:length(varargin)
    hyperface2elements(varargin{j},2) = -j;
end

idx = reshape( (1:nE)'*ones(1,nN),[],1) ~= ...
      hyperface2elements(element2hyperfaces(:),1);
element2neighbour = zeros(nE*nN,1);
element2neighbour(idx,1) = hyperface2elements(element2hyperfaces(idx),1);
element2neighbour(~idx,1) = hyperface2elements(element2hyperfaces(~idx),2);
element2neighbour = reshape(element2neighbour,[],nN);                          