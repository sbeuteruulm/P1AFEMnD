function M = stima(vertices)
%STIMA   Computes element stiffness matrix for simplex.
%   M = STIMA(X) computes element stiffness matrix for simplex.
%   The coordinates of the vertices are stored in X. 
%   X has dimension (d+1) x d. 
%
%   This routine should not be modified.

d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
M = abs(det([ones(1,d+1);vertices'])) * G * G' / factorial(d);
