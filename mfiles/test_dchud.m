% load test_dchud
N = 10;
F = randn(N, N);
A = F' * F / N;
x = randn(N, 1);
x = sign(x) .* x .^ 2;

% fprintf('Calculate new matrix:\n'); 
tic; At = A + x * x'; t2 = toc; 
% fprintf('elapsed time %g seconds.\n', t2);

% fprintf('Original cholesky:\n');
tic; L = chol(A); t2 = toc;
% fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
% norm(A - L' * L, 'fro'), norm(A - L' * L, 'fro') / norm(A, 'fro'), t2);

fprintf('Call chol:\n'); 
tic; Lt = chol(At); t2 = toc; 
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Lt' * Lt, 'fro'), norm(At - Lt' * Lt, 'fro') / norm(At, 'fro'), t2);

fprintf('Call cholupdate:\n'); 
tic; Ltm = cholupdate(L, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Ltm' * Ltm, 'fro'), norm(At - Ltm' * Ltm, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dcud up:\n'); 
tic; [Ltu, info] = dchud(L, 1, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltu' * Ltu, 'fro'), norm(At - Ltu' * Ltu, 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dcud lo:\n'); 
tic; [Ltl, info] = dchud(L', 0, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltl * Ltl', 'fro'), norm(At - Ltl * Ltl', 'fro') / norm(At, 'fro'), t2, info);
