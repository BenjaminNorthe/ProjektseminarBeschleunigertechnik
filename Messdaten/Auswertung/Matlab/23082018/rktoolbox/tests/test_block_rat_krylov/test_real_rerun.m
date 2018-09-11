function check = test_real_rerun
% This function tests rerunning when the real variant is used. This
% function also tests iterative refinement.

N = 80;
A = gallery('grcar', N, 3);
b = ones(N,2);
b(40,2) = 0;
xi = [0, 1-1i, 1+1i, 2+6i, 2-6i, inf];

param.real = 1;
param.refinement = 1;

tol = 1e-13;

[V, K, H] = rat_krylov(A, b, xi, param);
nrm1 = norm(A*V*K-V*H,'fro')/norm(H,'fro');
nrm2 = norm(V'*V-eye(size(V, 2)),'fro');

[Vr] = rat_krylov(A, b/norm(b), K ,H, param);
nrm3 = norm(A*Vr*K-Vr*H,'fro')/norm(H,'fro');

[Vr] = rat_krylov(A, rand(N,2), K ,H, param);
nrm4 = norm(A*Vr*K-Vr*H,'fro')/norm(H,'fro');

check = [nrm1 nrm2 nrm3 nrm4] < tol;
end