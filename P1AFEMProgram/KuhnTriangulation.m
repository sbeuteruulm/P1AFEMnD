function elements = KuhnTriangulation(n)
% splits a cube into simplices
% 
% elements = KuhnTriangulation(n) creates Kuhn simplices od a cube of
% dimension n. A cube is split into n! simplices that are defined by the
% vertices. The function creates a numbering of the vertices and defines
% the simplices by their vertex numbers. Array elements contains a Kuhn
% simplex in each row.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S.Beuter, S. Funken 18-10-22

A = perms(1:n);
pow2 = 2.^(0:n-1);
elements = cumsum([ones(size(A,1),1),pow2(A)],2);
