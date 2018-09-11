function check = test_complete_deflation()
% This code tests complete deflation when the whole space becomes
% invariant.

N = 20;
A = gallery('moler', N);
b = zeros(N,10);
b(1,1) = 1;
b(3,2) = 1;
b(5,3) = 1;
b(7,4) = 1;
b(9,5) = 1;
b(11,6) = 1;
b(13,7) = 1;
b(15,8) = 1;
b(17,9) = 1;
b(19,10) = 1;
xi = [inf, inf, inf];

tol = 1e-13;

s = warning; % turn of warning for lucky breakdown
warning('off')
[V, K, H] = rat_krylov(A, b, xi);
warning(s)

% This commented out norm is very large due to the deflated columns
% remaining in the decomposition so we maintain the Hessenberg pencil
% structure.
% nrm1 = norm(A*V*K-V*H, 'fro')/norm(H, 'fro')

nrm2 = norm(A*V(:,1:20)*K(1:20,1:10)-V(:,1:20)*H(1:20,1:10), 'fro')/norm(H(1:20,1:10), 'fro');
nrm3 = norm(V(:,1:20)'*V(:,1:20) - eye(20), 'fro');

check = [nrm2 nrm3] < tol;
end
