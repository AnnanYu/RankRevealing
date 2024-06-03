function [U,VT,A11inv,E22] = recompute(A,k)
    [m,n] = size(A);
    A11inv = inv(A(1:k,1:k));
    VT = A(1:k,1:k) \ A(1:k,(k+1):n); U = A((k+1):m,1:k) / A(1:k,1:k);
    E22 = A((k+1):m,(k+1):n) - A((k+1):m,1:k) * VT;
end