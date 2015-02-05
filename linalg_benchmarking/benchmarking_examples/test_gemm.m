%%
N = 200;
A = randn(5 * N, 10 * N); 
B = randn(10 * N, 50 * N);
At = A';
Bt = B';
numRuns = 100;

%%
% test standard
%
fprintf('Test product: A * B\n'); 
% fprintf('Matlab:'); 
% t2 = 0;
% for iter = 1:numRuns,
% 	tic; C = A * B; t2 = t2 + toc;
% end;
% fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL:\n');
tic; C1 = gemm_mkl_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma:\n');
tic; C2 = gemm_arma_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze:\n');
tic; C3 = gemm_blaze_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);

%%
% test left transpose
%
fprintf('Test product: A'' * B\n'); 
% fprintf('Matlab:'); 
% t2 = 0;
% for iter = 1:numRuns,
% 	tic; C = At' * B; t2 = t2 + toc;
% end;
% fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL:\n');
tic; C1 = gemm_mkl_test(At, B, 1, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma:\n');
tic; C2 = gemm_arma_test(At, B, 1, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze:\n');
tic; C3 = gemm_blaze_test(At, B, 1, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);

%%
% test right transpose
%
fprintf('Test product: A * B''\n'); 
% fprintf('Matlab:'); 
% t2 = 0;
% for iter = 1:numRuns,
% 	tic; C = A * Bt'; t2 = t2 + toc;
% end;
% fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL:\n');
tic; C1 = gemm_mkl_test(A, Bt, 0, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma:\n');
tic; C2 = gemm_arma_test(A, Bt, 0, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze:\n');
tic; C3 = gemm_blaze_test(A, Bt, 0, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);

%%
% test both transpose
%
fprintf('Test product: A'' * B''\n'); 
% fprintf('Matlab:'); 
% t2 = 0;
% for iter = 1:numRuns,
% 	tic; C = At' * Bt'; t2 = t2 + toc;
% end;
% fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL:\n');
tic; C1 = gemm_mkl_test(At, Bt, 1, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma:\n');
tic; C2 = gemm_arma_test(At, Bt, 1, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze:\n');
tic; C3 = gemm_blaze_test(At, Bt, 1, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);