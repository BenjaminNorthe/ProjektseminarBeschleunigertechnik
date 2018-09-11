function check = test_basic_arnoldi()

N = 100;
A = gallery('tridiag',N);
B = gallery('grcar', N, 3);
b = ones(N, 2);
b(1,1)   = 0;
b(4,1)   = -1;
b(N,1) = 0;
b(16,2) = -10;
b(17,2) = 3;

xi = [-3, 2, 3i, 8+12i, inf];

tol = 1e-13;

[V, K, H] = rat_krylov(A, b, xi);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V - eye(size(V,2)),'fro');

param.continuation = 'last';
[V, K, H] = rat_krylov(A, b, xi, param);
nrm3 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm4 = norm(V'*V - eye(size(V,2)),'fro');

param.continuation = 'ruhe';
param.refinement = 1;
[VV, KK, HH] = rat_krylov(A, B, b, xi, param);
nrm5 = norm(A*VV*KK - B*VV*HH, 'fro')/norm(HH, 'fro');
nrm6 = norm(VV'*VV - eye(size(VV,2)), 'fro');

check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6] <tol;
end