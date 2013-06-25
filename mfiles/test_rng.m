% load test_dchud
N = 10000;
M = 2 * N;
showVisualization = 0;

%%
% test normal
fprintf('Matlab generator:\n'); 
tic; R = randn(M, N); t2 = toc;
fprintf('elapsed time %g seconds.\n', t2);
if (showVisualization == 1)
    figure; hist(R(:), -6:0.1:6);
end;

fprintf('Call drnorm up:\n'); 
tic; Rd = drnorm(M, N); t2 = toc;
fprintf('elapsed time %g seconds.\n', t2);
if (showVisualization == 1)
    figure; hist(Rd(:), -6:0.1:6);
end;

%%
% test uniform
fprintf('Matlab generator:\n'); 
tic; U = rand(M, N); t2 = toc;
fprintf('elapsed time %g seconds.\n', t2);
if (showVisualization == 1)
    figure; hist(U(:), 0:0.1:1);
end;

fprintf('Call drunif up:\n'); 
tic; Ud = drunif(M, N); t2 = toc;
fprintf('elapsed time %g seconds.\n', t2);
if (showVisualization == 1)
    figure; hist(U(:), 0:0.1:1);
end;
