% Test script of RRQR and gammaQR
rng(37);
m = 100; n = 100; k = 50; gamma = 2;
A = rand(m,n);

% Compute gamma before any swap
[rat,~,~] = gammaQR(A,k);
fprintf('gamma with a uniformly random pivoting is %.3e\n',rat);

% Compute rank-revealing QR
[Q,R,Pi,Ak,swaps,t1,t2] = RRQR(A,k,gamma);

% Compute gamma after swaps
[rat,~,~] = gammaQR(A*Pi,k);
fprintf('gamma with a rank-revealing pivoting is %.3e\n',rat);