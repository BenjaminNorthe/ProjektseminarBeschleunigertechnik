function check = test_arnoldi_deflation()
% Tests deflation

% This first test is an arnoldi deflation, due to the choice of starting
% vector
N = 80;
A = gallery('tridiag',N);
b = ones(N,1);
b(40,1) = 0;

b = [b, (A^2)*b];
xi = [inf, inf, inf, inf, inf];
tol = 1e-13;

param.deflation_tol = 1e-10;
[V, K, H, out] = rat_krylov(A, b, xi, param);
nrm1 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm2 = norm(V'*V -eye(size(V, 2)), 'fro');
nrm3 = (out.blocksizes == [2,2,1,1,1,1]);

% Previously used param.continuation = 'ruhe', by default. Try one
% with param.continuation = 'last'.
param.continuation = 'last';
[V, K, H, outt] = rat_krylov(A, b, xi, param);
nrm4 = norm(A*V*K - V*H, 'fro')/norm(H, 'fro');
nrm5 = norm(V'*V -eye(size(V, 2)), 'fro');
nrm6 = (outt.blocksizes == [2,2,1,1,1,1]);

%%
% This next example uses a variant of the circulant matrix to shift away a
% non zero entry
A = sparse(N,N);
A(2:N+1:end) = 1;
A(1,N) = 0;

b = zeros(N,3);
b(1,1) = 1;
b(5,2) = 1; % after 4 iterations, this should deflate
b(N-2,3) = 1; % after 3 iterations, this should deflate

[VV, KK, HH, outtt] = rat_krylov(A, b, xi, param);

nrm7 = norm(A*VV*KK - VV*HH, 'fro')/norm(HH, 'fro');
nrm8 = norm(VV'*VV -eye(size(VV, 2)), 'fro');
nrm9 = (outtt.blocksizes == [3,3,3,2,1,1]);


check = [[nrm1 nrm2 nrm4 nrm5 nrm7 nrm8] < tol, nrm3 nrm6 nrm9];
end
