N = 10;
isUp = 0;
aSign = -1;

A = randn(N,N);
A = A'*A;
[U V] = eig(A);
A = U*sqrt(V)*U';
a = aSign * abs(randn(1,1)) / 10;
v = randn(N,1);
L = chol(A)';
if (isUp == 1), L = L'; end;
save temp A v L N a
L2 = chol(A + a * v * v')';
if (isUp == 1), L2 = L2'; end;
dbmex on
[L1 stat] = test_cholesky(L, isUp, v, a);
norm(L1-L2), stat
