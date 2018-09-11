function check = test_real_extend
% Test extending a  block rational Krylov space using the real variant.

N= 80;
A = gallery('grcar',N,3);
b = ones(N,2);
b(20,1) = 0;
b(60,2) = 0;

xi = [1+2i, 1-2i, 4+3i, 4-3i, -6+7i, -6-7i];

tol = 2e-13;
param.real = 1;

[V, K, H] = rat_krylov(A, b, xi, param);

nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
[VV, KK, HH] = rat_krylov(A, b, xi(1:4), param);

% Number of vectors required for extending
param.extend = 2;
[VVV, KKK, HHH] = rat_krylov(A, VV, KK, HH, xi(5:6),param);

% Extending a block Krylov decomposition uses 'CGS' as apposed to 'MGS' as
% previous block structure is unknown. This leads to slight column permutation
% differences in the rank revealing qr decomposition. 
nrm2 = norm(V*(V'*VVV) - VVV, 'fro')/norm(V, 'fro');
% nrm = norm(K - KKK, 'fro')/norm(K, 'fro');
% nrm = norm(H - HHH, 'fro')/norm(H, 'fro');
nrm3 = norm(A*VVV*KKK - VVV*HHH, 'fro')/norm(HHH, 'fro');

B = 0.5*gallery('moler', N) + eye(N);

[V, K, H] = rat_krylov(A, B, b, xi, param);
nrm4 = norm(A*V*K-B*V*H, 'fro')/norm(H,'fro');

[VV, KK, HH] = rat_krylov(A, B, b, xi(1:4), param);
[VVV, KKK, HHH] = rat_krylov(A, B, VV, KK, HH, xi(5:6), param);

nrm5 = norm(V*(V'*VVV) - VVV, 'fro')/norm(V, 'fro');
% nrm = norm(K - KKK, 'fro')/norm(K, 'fro');
% nrm = norm(H - HHH, 'fro')/norm(H, 'fro');
nrm6 = norm(A*VVV*KKK - B*VVV*HHH, 'fro')/norm(HHH, 'fro');

check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6] < tol;
end