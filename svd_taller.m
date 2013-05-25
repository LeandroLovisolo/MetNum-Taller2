function [U, S, V] = svd_taller(A)

S = diag(flipud(sqrt(eig(A'*A))));
[U LAMBDA] = eig(A*A');
[V LAMBDA] = eig(A'*A);
U = fliplr(U);
V = fliplr(V);