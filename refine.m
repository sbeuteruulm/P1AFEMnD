function s = refine(e)
% Performance of closure step of adaptive Newest Vertex Bisection
% refinement.
%
% s = refine(e) for a given element number e that shall be refined an array
% s with the numbers of elements that have to be refined additionally to
% keep the mesh conforming is returned.
%
% Comments:
%   Necessary mesh information is given as global arrays. For a marked
%   element any neighbour sharing the refinement edge has to be refined
%   as well. If two elements do not have the refinement edge in common,
%   additional refinements may be necessary. Then, bisection of allelements
%   around a refinement edge is done by the call of createSimultBisection.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S. Beuter, S. Funken 18-10-22

global E E2N nE 

nD = size(E,2) - 1; nN = size(E,2);
K = [];  F = e;
FK = false(nE,1);  FK(e) = true;    
while ~isempty(F)
  Fnew = [];
  for e1 = F(:)'
    e1a = E(e1,1); e1b = E(e1,nN);
    E2Nloc = E2N(e1,:);
    for j = 2 : nD
      e2 = E2Nloc(j);
      if e2 > 0 &&  ~FK(e2)
        %*** e2 has same refinement edge as e1
        e2a = E(e2,1); e2b = E(e2,nN);
        if  ((e1a == e2a) && (e1b == e2b)) || ...
            ((e1a == e2b) && (e1b == e2a))  
          Fnew = [Fnew,e2];
        else
          s = refine(e2);
          FK(nE) = false;            
          % add to Fnew the child of e2 that is a neighbour of e1                     
          if ismember(s(1),E2N(e1,:))
            Fnew = [Fnew,s(1)];
          else
            Fnew = [Fnew,s(2)];
          end
        end  
      end
    end
  end
  K = [K,F];
  F = Fnew;
  FK(F) = true; 
end
s = [K(1),nE+1];
K = unique(K);  i = find(K==e);
K = [K(i:end),K(1:i-1)];
createSimultBisection(K);