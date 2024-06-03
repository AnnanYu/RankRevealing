% Test script of RRLU and gammaLU
rng(37);
m = 100; n = 100; k = 50; gamma = 3;
A = rand(m,n);

% Compute gamma before any swap
[rat,~,~,~,~,~] = gammaLU(A,k);
fprintf('gamma with a uniformly random pivoting is %.3e\n',rat);

% Compute rank-revealing LU
[L,U,Pi1,Pi2,Ak,swaps,t1,t2] = RRLU(A,k,gamma);

% Compute gamma after swaps
[rat,~,~,~,~,~] = gammaLU(Pi1*A*Pi2,k);
fprintf('gamma with a rank-revealing pivoting is %.3e\n',rat);