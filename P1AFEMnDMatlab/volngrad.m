function varargout = volngrad(coordinates,elements)
% Calculation of volumes of each element and gradients of nodal interpolant
% functions
%
% varargout = volngrad(coordinates,elements) calculates the volume and
% gradient of the basis functions for a given simplex. The latter is an
% optional output variable. Accordingly, varargout is a cell array that
% contains either the volume only or the volume and an array for the
% gradients. The volumes are given as a vector and the gradients are cell 
% arrays themselves, one for each element. The input array coordinates
% saves for any mesh node the respective coordinates in a row. In the same
% way, the nodes of each element are given in the respective row of
% elements.
%
%Comments:
%   A Cholesky decomposition of the matrices containing the differences of
%   the element vertices is calculated. The advantage is that the volumes
%   can be calculated as a byproduct.
%
%Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
%Authors:
%   S. Beuter, S. Funken 18-10-22

nD = size(coordinates,2);
nN = size(elements,2);
nE = size(elements,1);
switch nN
  case 2
     varargout{1} = abs(coordinates(elements(:,1))-coordinates(elements(:,2)));
     if nargout == 2
         grad{1} = 1./(coordinates(elements(:,1))-coordinates(elements(:,2)));
         grad{2} = 1./(coordinates(elements(:,2))-coordinates(elements(:,1)));
         varargout{2} = grad;
     end
  otherwise    
    %*** Compute Gram matrix L
    idx = repmat(1:nN-1,nN-1,1);
    ii = idx(logical(tril(idx,0))); idx = idx';
    jj = idx(logical(tril(idx,0)));
    grad{nN,1} = coordinates(elements(:,nN),:);
    for k = 1:nN-1
      grad{k,1} = coordinates(elements(:,k),:) - grad{nN,1}; 
    end  
    L = zeros(nE,nN*(nN-1)/2);
    for k = 1:length(ii)
      L(:,k)=L(:,k)+dot(grad{ii(k)},grad{jj(k)},2);
    end
    %*** Compute Cholesky decomposition of L
    ptr = [1,1+cumsum(nN-1:-1:1)];
    L(:,ptr(1)) = sqrt(L(:,ptr(1)));
    for k = 1:nN-2
      L(:,ptr(k)+1:ptr(k+1)-1) =  L(:,ptr(k)+1:ptr(k+1)-1) ...
                             ./ L(:,ptr(k)*ones(ptr(k+1)-ptr(k)-1,1));
      for j = k+1:nN-1
        L(:,ptr(j):ptr(j+1)-1) = L(:,ptr(j):ptr(j+1)-1) ...
            - L(:,ptr(k)+j-k:ptr(k+1)-1) .* L(:,ones(1,nN-j)*(ptr(k)+j-k));
      end             
      L(:,ptr(k+1)) = sqrt(L(:,ptr(k+1)));
    end 
    %*** Compute volume
    varargout{1} = prod(L(:,ptr(1:nN-1)),2)/factorial(nN-1);
    if nargout == 2
      %*** Compute L y = A'
      for j = 1:nN-2
        grad{j} = grad{j}./L(:,ptr(j)*ones(1,nD));
        for k = j+1:nN-1
          grad{k} = grad{k} - grad{j}.*L(:,ones(1,nD)*(ptr(j)+k-j));  
        end
      end
      grad{nN-1}=grad{nN-1}./L(:,ptr(nN-1)*ones(1,nD));
      %*** Compute L' * x = y;
      grad{nN-1}=grad{nN-1}./L(:,ptr(nN-1)*ones(1,nD));
      grad{nN} = -grad{nN-1};
      for j = nN-2:-1:1
        for k = j+1:nN-1
          grad{j} = grad{j} - grad{k} .* L(:,ones(1,nD)*(ptr(j)+k-j));  
        end
        grad{j} = grad{j} ./ L(:,ptr(j*ones(1,nD)));
        grad{nN} = grad{nN} - grad{j};
      end
      varargout{2} = grad;
    end
end