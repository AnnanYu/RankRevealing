function [gamma,i,s] = gammaQR(A,k)
% gammaQR    Compute the local maximum volume ratio in a QR factorization
% 
% [gamma,i,s] = gammaQR(A,k) returns the value
%   of rat_lwb = max(vol(tilde(A_1)) / vol(A_1)), where A_1 = A(:,1:k)
%   and tilde(A_1) ranges over all m-by-k matrices
%   obtained from A_1 by at most a single column swap.
%   In addition, the bound can be achieved
%   by swapping column i with column s.
    
    gamma = 1; i = 1; s = 1;
    [~,R] = qr(A); i = 0; s = 0;
    R11 = R(1:k,1:k); R12 = R(1:k,(k+1):end); R22 = R((k+1):end,(k+1):end);
    Nu = R11 \ R12; A1 = inv(R11'*R11); A2 = R22'*R22;
    for j = 1:k
        for t = (k+1):size(A,2)
            ratio = sqrt(Nu(j,t-k)^2 + A1(j,j) * A2(t-k,t-k));
            if ratio > gamma
                gamma = ratio; i = j; s = t;
            end
        end
    end
end