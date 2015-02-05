%%
N = 200;
A = randn(10 * N, 10 * N); 
A = (A + A') / 2;
B = randn(10 * N, 10 * N);
numRuns = 100;

%%
% test left
%
fprintf('Test product: A * B\n'); 
% fprintf('Matlab:'); 
% t2 = 0;
% for iter = 1:numRuns,
% 	tic; C = A * B; t2 = t2 + toc;
% end;
% fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL up:\n');
tic; C1 = symm_mkl_test(A, B, 1, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL lo:\n');
tic; C2 = symm_mkl_test(A, B, 1, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL gemm:\n');
tic; C3 = gemm_mkl_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma up:\n');
tic; C4 = symm_arma_test(A, B, 1, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma lo:\n');
tic; C5 = symm_arma_test(A, B, 1, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma gemm:\n');
tic; C6 = gemm_arma_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze:\n');
tic; C7 = symm_blaze_test(A, B, 1, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze gemm:\n');
tic; C8 = gemm_blaze_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);

%%
% test right
%
fprintf('Test product: B * A\n'); 
% fprintf('Matlab:'); 
% t2 = 0;
% for iter = 1:numRuns,
% 	tic; C = B * A; t2 = t2 + toc;
% end;
% fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL up:\n');
tic; C1 = symm_mkl_test(A, B, 0, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL lo:\n');
tic; C2 = symm_mkl_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('MKL gemm:\n');
tic; C3 = gemm_mkl_test(B, A, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma up:\n');
tic; C4 = symm_arma_test(A, B, 0, 1, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma lo:\n');
tic; C5 = symm_arma_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Arma gemm:\n');
tic; C6 = gemm_arma_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze:\n');
tic; C7 = symm_blaze_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);
fprintf('Blaze gemm:\n');
tic; C8 = gemm_blaze_test(A, B, 0, 0, numRuns); t2 = toc;
fprintf('\telapsed time %g seconds.\n', t2 / numRuns);