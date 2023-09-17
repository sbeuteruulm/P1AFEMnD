function [hyperface2nodes,element2hyperfaces,varargout] = ...
                                    provideHyperFaceData(elements,varargin)
nE = size(elements,1);
nN = size(elements,2);
nB = length(varargin);
nBF = zeros(1,nB);
for j = 1:nB
    nBF(j) = size(varargin{j},1);
end
ptr = cumsum([nE*nN,nBF]);
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
    varargout{j} = reshape(J(ptr(j)+1:ptr(j+1)),nBF(j),1);
end