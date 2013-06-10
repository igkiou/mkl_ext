% load test_dchud
N = 1000;
K = 50;
F = randn(N, N);
A = F' * F / N;
Y = randn(N, K);
Y = sign(Y) .* Y .^ 2;
YN = Y';
evec = eig(A);
apos = (abs(min(evec)) / 2) ^ 2;
aneg = - apos;
X = Y * sqrt(apos);
XN = X';

tic; L = chol(A); t2 = toc;

%
% positive tests first
%

% fprintf('Calculate new matrix:\n'); 
tic; At = A + X * X'; t2 = toc; 
% fprintf('elapsed time %g seconds.\n', t2);

fprintf('Call chol:\n'); 
tic; Lt = chol(At); t2 = toc; 
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(At - Lt' * Lt, 'fro'), norm(At - Lt' * Lt, 'fro') / norm(At, 'fro'), t2);

fprintf('Call dchrk up:\n'); 
tic; [Ltru, info] = dchrk(L, 1, Y, 1, apos); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltru' * Ltru, 'fro'), norm(At - Ltru' * Ltru, 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dchrk lo:\n'); 
Ltran = L';
tic; [Ltrl, info] = dchrk(Ltran, 0, Y, 1, apos); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltrl * Ltrl', 'fro'), norm(At - Ltrl * Ltrl', 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dchrk up:\n'); 
tic; [Ltru, info] = dchrk(L, 1, YN, 0, apos); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltru' * Ltru, 'fro'), norm(At - Ltru' * Ltru, 'fro') / norm(At, 'fro'), t2, info);

fprintf('Call dchrk lo:\n'); 
Ltran = L';
tic; [Ltrl, info] = dchrk(Ltran, 0, YN, 0, apos); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds,\tstatus %d.\n',...
norm(At - Ltrl * Ltrl', 'fro'), norm(At - Ltrl * Ltrl', 'fro') / norm(At, 'fro'), t2, info);
