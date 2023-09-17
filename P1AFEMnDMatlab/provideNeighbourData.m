function [element2neighbour,hyperface2elements] = ...
                          provideNeighbourData(element2hyperfaces,varargin)
% Identification of neighbouring elements, i.e. elements that share a
% hyperface.
%
% [element2neighbour,hyperface2elements] = ...
% provideNeighbourData(element2hyperfaces,varargin) provides information on
% neighbour relations. As input two arguments are expected: an array
% containing the elements given by their hyperface numbers saved row wise
% and a cell array varargin that gives the hyperfaces of all boundary
% types. For each element, the neighbour numbers are returned in the rows
% of array element2neighbour. In hyperface2elements each row contains the
% pair of elements that share a hyperface.
%
%Comments:
%   The neighbour relations are crucial for the realization of the adaptive
%   FEM. They are provided in this function and can be updated when a mesh
%   is refined.
%
%Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
%Authors:
%   S. Beuter, S. Funken 18-10-22

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