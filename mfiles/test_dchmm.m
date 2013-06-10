% load test_dchud
N = 10000;
M = 2 * N;
F = randn(N, N);
AN = F' * F / N;
F = randn(M, M);
AM = F' * F / M;
B = randn(M, N);
alpha = randn(1);
LM = chol(AM); 
LN = chol(AN); 

%
% test side = 'L'
%
fprintf('Standard product:\n'); 
tic; R = alpha * (AM * B); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(R - R, 'fro'), norm(R - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Matlab cholesky no parenthesis:\n'); 
tic; Rn = alpha * LM' * LM * B; t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rn - R, 'fro'), norm(Rn - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Matlab cholesky wrong parenthesis:\n'); 
tic; Rw = alpha * (LM' * LM) * B; t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rw - R, 'fro'), norm(Rw - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Matlab cholesky parenthesis:\n'); 
tic; Rm = alpha * (LM' * (LM * B)); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rm - R, 'fro'), norm(Rm - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Call dchmm up:\n'); 
tic; Ru = dchmm(LM, 1, B, 1, alpha); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Ru - R, 'fro'), norm(Ru - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Call dchmm lo:\n'); 
Lt = LM';
tic; Rl = dchmm(Lt, 0, B, 1, alpha); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rl - R, 'fro'), norm(Rl - R, 'fro') / norm(R, 'fro'), t2);

%
% test side = 'R'
%
fprintf('Standard product:\n'); 
tic; R = alpha * (B * AN); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(R - R, 'fro'), norm(R - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Matlab cholesky no parenthesis:\n'); 
tic; Rn = alpha * B * LN' * LN ; t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rn - R, 'fro'), norm(Rn - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Matlab cholesky wrong parenthesis:\n'); 
tic; Rw = alpha * B * (LN' * LN); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rw - R, 'fro'), norm(Rw - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Matlab cholesky parenthesis:\n'); 
tic; Rm = alpha * ((B * LN') * LN); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rm - R, 'fro'), norm(Rm - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Call dchmm up:\n'); 
tic; Ru = dchmm(LN, 1, B, 0, alpha); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Ru - R, 'fro'), norm(Ru - R, 'fro') / norm(R, 'fro'), t2);

fprintf('Call dchmm lo:\n'); 
Lt = LN';
tic; Rl = dchmm(Lt, 0, B, 0, alpha); t2 = toc;
fprintf('error %g,\trelative error %g,\telapsed time %g seconds.\n',...
norm(Rl - R, 'fro'), norm(Rl - R, 'fro') / norm(R, 'fro'), t2);
