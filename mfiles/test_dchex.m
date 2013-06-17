% load test_dchud
N = 10;
F = randn(N, N);
A = F' * F / N;
k = ceil(N * 0.2);
l = ceil(N * 0.7);
E = eye(N);
E(:, k:l) = E(:, [l k:(l - 1)]);

%
fprintf('Test job = 1.\n');
% fprintf('Calculate new matrix:\n'); 
tic; At = E' * A * E; t2 = toc; 
% fprintf('elapsed time %g seconds.\n', t2);

% fprintf('Original cholesky:\n');
tic; L = chol(A); t2 = toc;
% fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
% norm(A - L' * L, 'fro'), norm(A - L' * L, 'fro') / norm(A, 'fro'), t2);

fprintf('Call chol:\n'); 
tic; Lt = chol(At); t2 = toc; 
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Lt' * Lt, 'fro'), norm(At - Lt' * Lt, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dchex:\n'); 
tic; Ltm = dchex_mat(L, N, N,  k, l, [], 0, 0, 1); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Ltm' * Ltm, 'fro'), norm(At - Ltm' * Ltm, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dchex:\n'); 
tic; Ltu = dchex(L, k, l, 1); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Ltu' * Ltu, 'fro'), norm(At - Ltu' * Ltu, 'fro') / norm(At, 'fro'), t2);

%
fprintf('Test job = 2.\n');
% fprintf('Calculate new matrix:\n'); 
tic; At = E * A * E'; t2 = toc; 
% fprintf('elapsed time %g seconds.\n', t2);

% fprintf('Original cholesky:\n');
tic; L = chol(A); t2 = toc;
% fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
% norm(A - L' * L, 'fro'), norm(A - L' * L, 'fro') / norm(A, 'fro'), t2);

fprintf('Call chol:\n'); 
tic; Lt = chol(At); t2 = toc; 
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Lt' * Lt, 'fro'), norm(At - Lt' * Lt, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dchex:\n'); 
tic; Ltm = dchex_mat(L, N, N,  k, l, [], 0, 0, 2); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Ltm' * Ltm, 'fro'), norm(At - Ltm' * Ltm, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dchex:\n'); 
tic; Ltu = dchex(L, k, l, 2); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Ltu' * Ltu, 'fro'), norm(At - Ltu' * Ltu, 'fro') / norm(At, 'fro'), t2);
