% load test_dchud
N = 10000;
F = randn(N, N);
A = F' * F / N;
x = randn(N, 1);

% fprintf('Original cholesky:\n');
tic; L = chol(A); t2 = toc;
% fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
% norm(A - L' * L, 'fro'), norm(A - L' * L, 'fro') / norm(A, 'fro'), t2);

fprintf('Standard product:\n'); 
tic; y = A * x; t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(y - y, 'fro'), norm(y - y, 'fro') / norm(y, 'fro'), t2);

fprintf('Matlab cholesky no parenthesis:\n'); 
tic; yn = L' * L * x; t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(yn - y, 'fro'), norm(yn - y, 'fro') / norm(y, 'fro'), t2);

fprintf('Matlab cholesky wrong parenthesis:\n'); 
tic; yw = (L' * L) * x; t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(yw - y, 'fro'), norm(yw - y, 'fro') / norm(y, 'fro'), t2);

fprintf('Matlab cholesky parenthesis:\n'); 
tic; ym = (L' * (L * x)); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(ym - y, 'fro'), norm(ym - y, 'fro') / norm(y, 'fro'), t2);

fprintf('Call dchmv up:\n'); 
tic; yu = dchmv(L, 1, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(yu - y, 'fro'), norm(yu - y, 'fro') / norm(y, 'fro'), t2);

fprintf('Call dchmv lo:\n'); 
Lt = L';
tic; yl = dchmv(Lt, 0, x); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(yl - y, 'fro'), norm(yl - y, 'fro') / norm(y, 'fro'), t2);
