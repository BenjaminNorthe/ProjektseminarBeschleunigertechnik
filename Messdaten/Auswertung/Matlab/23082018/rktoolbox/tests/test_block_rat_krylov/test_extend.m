function check = test_extend
% Test extending a block rational Krylov space.

A = rand(50,50) + 1i*rand(50,50);
b = rand(50,3) + 1i*rand(50,3);
xi = [5, 7, 4+3i, 2, inf, 0, 9];

tol = 1e-13;

[V, K, H] = rat_krylov(A, b, xi);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
[VV, KK, HH] = rat_krylov(A, b, xi(1:5));

% Number of vectors required for extending
param.extend = 3;
param.orth = 'CGS';
[VVV1, KKK1, HHH1] = rat_krylov(A, VV, KK, HH, xi(6:7), param);

% Extending a block Krylov decomposition uses 'CGS' and 'MGS' differently as
% previous block structure is unknown. This leads to slight column permutation
% differences in the rank revealing qr decomposition. 
P = V'*VVV1;
nrm2 = norm(V*P - VVV1, 'fro')/norm(V, 'fro');
nrm3 = norm(A*VVV1*KKK1 - VVV1*HHH1, 'fro')/norm(HHH1, 'fro');
% nrm = norm(K - KKK1, 'fro')/norm(K, 'fro');
% nrm = norm(H - HHH1, 'fro')/norm(H, 'fro');

param.orth = 'MGS';
[VVV2, KKK2, HHH2] = rat_krylov(A, VV, KK, HH, xi(6:7), param);

P = V'*VVV2;
nrm4 = norm(V*P - VVV2, 'fro')/norm(V, 'fro');
nrm5 = norm(A*VVV2*KKK2 - VVV2*HHH2, 'fro')/norm(HHH2, 'fro');
% nrm = norm(K - KKK, 'fro')/norm(K, 'fro');
% nrm = norm(H - HHH, 'fro')/norm(H, 'fro');

B = rand(50,50) + 1i*rand(50,50);

[V, K, H] = rat_krylov(A, B, b, xi);
nrm6 = norm(A*V*K-B*V*H, 'fro')/norm(H,'fro');
[VV, KK, HH] = rat_krylov(A, B, b, xi(1:5));

[VVV, KKK, HHH] = rat_krylov(A, B, VV, KK, HH, xi(6:7), param);

nrm7 = norm(V*(V'*VVV) - VVV, 'fro')/norm(V, 'fro');
nrm8 = norm(A*VVV*KKK - B*VVV*HHH, 'fro')/norm(HHH, 'fro');
% nrm = norm(K - KKK, 'fro')/norm(K, 'fro');
% nrm = norm(H - HHH, 'fro')/norm(H, 'fro');

check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7 nrm8] < tol;
end