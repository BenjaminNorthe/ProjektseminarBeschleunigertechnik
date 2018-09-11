function check = test_rerun
% Test rerunning a block rational Arnoldi decomposition.

N = 70;
A = gallery('tridiag', N);
b = ones(N,3);
b(15,2) = 0;
b(30,3) = 5;
xi = [0, 1, 2, 3, inf];

tol = 1e-13;

[V, K, H] = rat_krylov(A, b, xi);

nrm1 = norm(A*V*K-V*H,'fro');
nrm1 = nrm1/(norm(A,'fro')*norm(V,'fro')*norm(K,'fro')+norm(V,'fro')*norm(H,'fro'));
nrm2 = norm(V'*V-eye(size(V, 2)),'fro');

Vr = rat_krylov(A, b/norm(b), K, H);

nrm3 = norm(V*(V'*Vr)-Vr)/length(xi);
nrm4 = norm(A*Vr*K-Vr*H,'fro');
nrm4 = nrm4/(norm(A,'fro')*norm(V,'fro')*norm(K,'fro')+norm(V,'fro')*norm(H,'fro'));

c = rand(N,3);

Vrr = rat_krylov(A, c/norm(c), K, H);

nrm5 = norm(A*Vrr*K-Vrr*H,'fro');
nrm5 = nrm4/(norm(A,'fro')*norm(V,'fro')*norm(K,'fro')+norm(V,'fro')*norm(H,'fro'));

check = [nrm1 nrm2 nrm3 nrm4 nrm5] < tol;
end