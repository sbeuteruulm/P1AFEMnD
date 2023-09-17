function A = convertdec(n,b)
% converts integers 0,1,2,...,b^n-1 from decimal
% to number system with base b (reverse order)
A = zeros(b^n,n);
idx = (0:b^n-1)';
%for i=1:n
for i = n:-1:1
  A(:,i) = rem(idx,b);
  idx = floor(idx/b);
end