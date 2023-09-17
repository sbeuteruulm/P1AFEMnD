function elements = KuhnTriangulation(n)
A = perms(1:n);
pow2 = 2.^(0:n-1);
elements = cumsum([ones(size(A,1),1),pow2(A)],2);
