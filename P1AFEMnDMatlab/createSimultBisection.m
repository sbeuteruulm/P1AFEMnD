function createSimultBisection(varargin)
% Simultaneous bisection of a marked element and its neighbours sharing
% the refinement edge.
% 
% createSimultBisection(varargin) updates a globally defined mesh by
% bisection of a marked element and its neighbours sharing a refinement
% edge. An element number of the simplex that shall be refined is expected
% as input varargin.
% 
% Comments:
%   All elements around the refinement edge can be bisected
%   simultaneously. The element array gets enlarged. Since the
%   neighbour information is crucial for the refinement process, it has to
%   be updated for all bisected elements, too. With each call of
%   createSimultBisection one new mesh node is calculated and added to the
%   array coordinates.
% 
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
% 
% Authors:
%   S. Beuter, S. Funken 18-10-22

global C E L E2N i2fi nC nE                               
marked = varargin{end};
nR = length(marked);      % number of simplices to be refined
nN = size(E,2);           % number of nodes per simplex
% create new node
ae = E(marked(1),[1,end]);
C(nC+1,:) = sum(C(ae,:))/2; nC = nC+1;
% create new elements and update element2neighbour
elem = E(marked,:);
low = [elem(:,1  ),ones(nR,1)*(nC), elem(:,2:end-1)];
up  = [elem(:,end),ones(nR,1)*(nC), elem(:,2:end-1)];
neigh = E2N(marked,:);
n_low  = [nE+(1:nR)', neigh(:,end), neigh(:,2:end-1)];
i2fi(marked,1) = 1:nR;
tmp = neigh(:,2:end-1);
idx = find(tmp > 0);
tmp(idx) = nE+i2fi(tmp(idx));
n_up = [marked(:), neigh(:,1), tmp];
ind = find(elem(:,1) > elem(:,end));
lev = mod(L(marked),nN-1);
tmp=n_up(ind,2:end);n_up(ind,2:end)=n_low(ind,2:end);n_low(ind,2:end)=tmp;
% reestablish correct node order
for k = 0:nN-3 %nD-1
  idx = find( lev == k );
  up(idx,3+k:end) = up(idx,end:-1:3+k);
  n_up(idx,3+k:end) = n_up(idx,end:-1:3+k);
end
% update global arrays for elements and nieghours
tmp=n_up(ind,3:end);n_up(ind,3:end)=n_low(ind,3:end);n_low(ind,3:end)=tmp;
E2N(marked,:) = n_low;
E2N(nE+(1:nR)',:) = n_up;
tmp = up(ind,:); up(ind,:) = low(ind,:); low(ind,:) = tmp;

E(marked,:) = low;
E(nE+(1:nR)',:) = up;
% update element2neighbour for elements opposite the first node
ii = find(neigh(:,1)>0);
imat=E2N(neigh(ii,1),:)==reshape(marked(ii),[],1)*ones(1,nN);
[i,j,k] = find(imat);
tmp = n_low(:,1); tmp(ind) = n_up(ind,1);
for k=1:length(i)
  E2N(neigh(ii(i(k)),1),j(k)) = tmp(ii(i(k)));
end
% update element2neighbour for elements opposite the last node
ii = find(neigh(:,end)>0);
imat=E2N(neigh(ii,end),:)==reshape(marked(ii),[],1)*ones(1,nN);
[i,j,k] = find(imat);
tmp = n_up(:,1); tmp(ind) = n_low(ind,1);
for k=1:length(i)
  E2N(neigh(ii(i(k)),end),j(k)) = tmp(ii(i(k)));
end
% update level
L([marked,nE+(1:nR)]) = [L(marked),L(marked)] + 1;
nE = nE + nR;

