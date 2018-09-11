function check = test_extend_deflation
% Test extending a block rational Krylov space when deflation has occurred.

N = 80;
A = gallery('tridiag',N);
b = ones(N,1);
b(floor(N/2),1) = 0;

b = [b, (A^3)*b];
xi = [inf, inf, inf, inf, inf, inf, inf];

param.deflation_tol = 1e-10;

tol = 1e-13;

[V, K, H, out] = rat_krylov(A, b, xi, param);

% Fat decomposition
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
% Slim decomposition
nrm2 = norm(A*V*K(:, out.column_deflation) - V*H(:,out.column_deflation), 'fro')/norm(H, 'fro');

[VV, KK, HH, outt] = rat_krylov(A, b, xi(1:5), param);

% Number of vectors required for extending
param.extend = 2;

% Extending a block Krylov decomposition uses 'CGS' as apposed to 'MGS' as
% previous block structure is unknown. This leads to slight column permutation
% differences in the rank revealing qr decomposition. 

% In delfation case, must extend with slim decomposition, or change continutation vector,
% otherwise 'ruhe' continuation vector will not work.
[VVV, KKK, HHH, outtt] = rat_krylov(A, VV, KK(:, outt.column_deflation), HH(:, outt.column_deflation), xi(6:7), param);

P = (V'*VVV);
nrm3 = norm(V*P - VVV, 'fro')/norm(V, 'fro');
% nrm = norm(K(:,out.column_deflation) - KKK(:,outtt.column_deflation), 'fro')/norm(K(:, out.column_deflation), 'fro');
% nrm = norm(H(:,out.column_deflation) - HHH(:,outtt.column_deflation), 'fro')/norm(H(:, out.column_deflation), 'fro');

% Fat decomposition
nrm4 = norm(A*VVV*KKK - VVV*HHH, 'fro')/norm(HHH, 'fro');
% Slim decomposition
nrm5 = norm(A*VVV*KKK(:, outtt.column_deflation) - VVV*HHH(:,outtt.column_deflation), 'fro')/norm(HHH, 'fro');

%%

B = gallery('grcar', N, 3);

[V, K, H, out] = rat_krylov(A, B, b, xi);
% Fat decomposition
nrm6 = norm(A*V*K - B*V*H, 'fro')/norm(H, 'fro');
% Slim decomposition
nrm7 = norm(A*V*K(:, out.column_deflation) - B*V*H(:,out.column_deflation), 'fro')/norm(H, 'fro');

[VV, KK, HH] = rat_krylov(A, B, b, xi(1:5));

[VVV, KKK, HHH] = rat_krylov(A, B, VV, KK, HH, xi(6:7), param);

P = (V'*VVV);
nrm8 = norm(V*P - VVV, 'fro')/norm(V, 'fro');
% nrm = norm(K - KKK, 'fro')/norm(K, 'fro');
% nrm = norm(H - HHH, 'fro')/norm(H, 'fro');
% Fat decomposition
nrm9 = norm(A*VVV*KKK - B*VVV*HHH, 'fro')/norm(HHH, 'fro');
% Slim decomposition
nrm10 = norm(A*VVV*KKK(:, outtt.column_deflation) - B*VVV*HHH(:,outtt.column_deflation), 'fro')/norm(HHH, 'fro');

check = [nrm1 nrm2 nrm3 nrm4 nrm5 nrm6 nrm7 nrm8 nrm9 nrm10] < tol;
end