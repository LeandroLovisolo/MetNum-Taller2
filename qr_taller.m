function [Q R] = qr_taller(S, V)

SV = S * V';

% Givens en (2,1)
x1 = SV(1,1);
x2 = SV(2,1);
n = norm([x1, x2]);
G1 = [x1/n x2/n 0; -x2/n x1/n 0; 0 0 1];
SV = G1 * SV;

% Givens en (3,1)
x1 = SV(1,1);
x2 = SV(3,1);
n = norm([x1, x2]);
G2 = [x1/n 0 x2/n; 0 1 0; -x2/n 0 x1/n];
SV = G2 * SV;

% Householder en (3,2)
x = SV(2:3,2);
y = [norm(x); 0];
u = x - y;
H = eye(2) - 2*u*u' / norm(u)^2;
H3 = [1 0 0; 0 H(1,1) H(1,2); 0 H(2,1) H(2,2)];
SV = H3 * SV;

Q = G1' * G2' * H3';
%Q = G1^-1 * G2^-1 * H3^-1;
R = SV;