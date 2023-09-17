function A = convertdec(n,b)
% converts integers 0,1,2,...,b^n-1 from decimal
% to number system with base b (reverse order)
%
% A = convertdec(n,b) converts all natural numbers 0 to n to a number
% system with base b. A is a matrix that contains the numbers 0 to n in
% each row in ascending order in the converted number system.
%
% Remark:
%   This program is a supplement to the paper 
%   >> Efficient P1-FEM for any space dimension in Matlab <<
%   by S. Beuter, and S. Funken. The reader should 
%   consult that paper for more information.   
%
% Authors:
%   S. Beuter, S. Funken 18-10-22

A = zeros(b^n,n);
idx = (0:b^n-1)';
%for i=1:n
for i = n:-1:1
  A(:,i) = rem(idx,b);
  idx = floor(idx/b);
end