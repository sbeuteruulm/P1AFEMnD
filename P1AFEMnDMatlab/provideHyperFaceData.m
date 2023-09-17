function [hyperface2nodes,element2hyperfaces,varargout] = ...
                                    provideHyperFaceData(elements,varargin)
% Calculation of additional mesh data. Especially, the hyperface data is
% collected.
%
% [hyperface2nodes,element2hyperfaces,varargout] = ...
% provideHyperFaceData(elements,varargin) Creation of a numbering for
% hyperfaces of a mesh described by their node numbers, which are given for
% each element in elements. Each hyperface is given by its nodes row wise
% in hyerface2elements. The hyperfaces of the simplices are assembled in
% element2hyperfaces where each row contains the hyperfaces of an element.
% According to the cell array varargin which contains all nodes of each
% boundary type, a cell array varargout is created that contains all
% hyperfaces of each boundary type.
%
% Comments:
%   For any boundary type the hyperfaces being part of it are also given in
%   an array.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S. Beuter, S. Funken 18-10-22

nE = size(elements,1);
nN = size(elements,2);
nB = length(varargin);
nBE = zeros(1,nB);
for j = 1:nB
    nBE(j) = size(varargin{j},1);
end
ptr = cumsum([nE*nN,nBE]);
% i'th hyperface is opposite the i'th node 
ordering = flipud(nchoosek(1:nN,nN-1));
%*** ascending sorting of nodes of all facets 
facets = sort(reshape(elements(:,ordering),[],nN-1),2);
%*** add boundary facets 
for j = 1:nB
    boundary = varargin{j};
    if ~isempty(boundary)
        facets(ptr(j)+1:ptr(j+1),:) = sort(boundary,2);
    end
end
%*** create unique facet numbering, element2hyperfaces, and boundary2hyperfaces
[hyperface2nodes,~,J] = unique(facets,'rows');
element2hyperfaces = reshape(J(1:nN*nE),[],nN);
varargout = cell(nB,1);
for j = 1:nB
    varargout{j} = reshape(J(ptr(j)+1:ptr(j+1)),nBE(j),1);
end