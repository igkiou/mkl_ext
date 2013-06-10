% load test_dchud
N = 10000;
F = randn(N, N);
A = F' * F / N;
y = randn(N, 1);
y = sign(y) .* y .^ 2;
evec = eig(A);
apos = (abs(min(evec)) / 2) ^ 2;
aneg = - apos;
x = y * sqrt(apos);

tic; L = chol(A); t2 = toc;

%
% positive tests first
%

% fprintf('Calculate new matrix:\n'); 
tic; At = A + x * x'; t2 = toc; 
% fprintf('elapsed time %g seconds.\n', t2);

fprintf('Call chol:\n'); 
tic; Lt = chol(At); t2 = toc; 
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Lt' * Lt, 'fro'), norm(At - Lt' * Lt, 'fro') / norm(At, 'fro'), t2);

fprintf('Call cholupdate:\n'); 
tic; Ltm = cholupdate(L, x, '+'); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Ltm' * Ltm, 'fro'), norm(At - Ltm' * Ltm, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dchud up:\n'); 
tic; [Ltu, info] = dchud(L, 1, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltu' * Ltu, 'fro'), norm(At - Ltu' * Ltu, 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dchud lo:\n'); 
Ltran = L';
tic; [Ltl, info] = dchud(Ltran, 0, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltl * Ltl', 'fro'), norm(At - Ltl * Ltl', 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dchr up:\n'); 
tic; [Ltru, info] = dchr(L, 1, y, apos); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltru' * Ltru, 'fro'), norm(At - Ltru' * Ltru, 'fro') / norm(At, 'fro'), t2, info);
fprintf('difference %g,\trelative difference %g\n',...
norm(Ltu - Ltru, 'fro'), norm(Ltu - Ltru, 'fro') / norm(Ltu, 'fro'));

fprintf('Call dchr lo:\n'); 
Ltran = L';
tic; [Ltrl, info] = dchr(Ltran, 0, y, apos); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltrl * Ltrl', 'fro'), norm(At - Ltrl * Ltrl', 'fro') / norm(At, 'fro'), t2, info);
fprintf('difference %g,\trelative difference %g\n',...
norm(Ltl - Ltrl, 'fro'), norm(Ltl - Ltrl, 'fro') / norm(Ltl, 'fro'));

%
% negative tests second
%

% fprintf('Calculate new matrix:\n'); 
tic; At = A - x * x'; t2 = toc; 
% fprintf('elapsed time %g seconds.\n', t2);

fprintf('Call chol:\n'); 
tic; Lt = chol(At); t2 = toc; 
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Lt' * Lt, 'fro'), norm(At - Lt' * Lt, 'fro') / norm(At, 'fro'), t2);

fprintf('Call cholupdate:\n'); 
tic; Ltm = cholupdate(L, x, '-'); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Ltm' * Ltm, 'fro'), norm(At - Ltm' * Ltm, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dchud up:\n'); 
tic; [Ltu, info] = dchdd(L, 1, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltu' * Ltu, 'fro'), norm(At - Ltu' * Ltu, 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dchud lo:\n'); 
Ltran = L';
tic; [Ltl, info] = dchdd(Ltran, 0, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltl * Ltl', 'fro'), norm(At - Ltl * Ltl', 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dchr up:\n'); 
tic; [Ltru, info] = dchr(L, 1, y, aneg); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltru' * Ltru, 'fro'), norm(At - Ltru' * Ltru, 'fro') / norm(At, 'fro'), t2, info);
fprintf('difference %g,\trelative difference %g\n',...
norm(Ltu - Ltru, 'fro'), norm(Ltu - Ltru, 'fro') / norm(Ltu, 'fro'));

fprintf('Call dchr lo:\n'); 
Ltran = L';
tic; [Ltrl, info] = dchr(Ltran, 0, y, aneg); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltrl * Ltrl', 'fro'), norm(At - Ltrl * Ltrl', 'fro') / norm(At, 'fro'), t2, info);
fprintf('difference %g,\trelative difference %g\n',...
norm(Ltl - Ltrl, 'fro'), norm(Ltl - Ltrl, 'fro') / norm(Ltl, 'fro'));
