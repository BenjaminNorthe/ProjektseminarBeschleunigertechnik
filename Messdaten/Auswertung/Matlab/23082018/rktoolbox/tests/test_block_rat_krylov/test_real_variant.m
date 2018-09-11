function check = test_real_variant()
% This function tests using the real variant for complex conjugate poles.

N = 60;
A = gallery('dorr', N);
A = A + ones(N); % A is fairly unconditioned.
b = ones(N,3);
b(4,1) = 0;
b(16,2) = 12;
b(47,3) = -4;
xi = [2, 1+1i, 1-1i, 4i, -4i, 0, inf];

tol = 1e-13;

param.real = 1;

[V,K,H, out] = rat_krylov(A, b, xi, param);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)), 'fro');

param.real = 0;

[VV,KK,HH, out] = rat_krylov(A, b, xi, param);
nrm3 = norm(A*VV*KK - VV*HH, 'fro')/norm(HH, 'fro');
nrm4 = norm(VV'*VV - eye(size(VV,2)), 'fro');

check = [nrm1 nrm2 nrm3 nrm4] < tol;
end
