function [Q,R,Pi,gamma,R12,A2] = CPQR(A,kf)
    X = A; [m,n] = size(A); Pi = eye(n);
    % if kf >= min(m,n)
    %     [Q, R, ~] = qr(A,0); gamma = 0; R12 = 0; A2 = 0; return
    % end
    gamma = vecnorm(A); [~,j] = max(gamma);
    Pi(:,[1,j]) = Pi(:,[j,1]); A(:,[1,j]) = A(:,[j,1]); gamma([1,j]) = gamma([j,1]);
    R = norm(A(:,1)); Q = A(:,1) / R;
    B = Q' * A(:,2:n); A(:,2:n) = A(:,2:n) - Q * B;
    gamma(1) = 0; gamma(2:n) = sqrt(gamma(2:n).^2 - B.^2);
    for k = 2:kf
        [~,j] = max(gamma);
        Pi(:,[k,j]) = Pi(:,[j,k]); A(:,[k,j]) = A(:,[j,k]); gamma([k,j]) = gamma([j,k]);
        B(:,[1,j-k+1]) = B(:,[j-k+1,1]);
        r = norm(A(:,k)); Q = [Q,A(:,k)/r]; R = [R,B(:,1);zeros(1,k-1),r];
        b = Q(:,k)' * A(:,(k+1):n); B(:,1) = []; B = [B;b];
        A(:,(k+1):n) = A(:,(k+1):n) - Q(:,k) * b;
        gamma(k) = 0; gamma((k+1):n) = sqrt(gamma((k+1):n).^2 - b.^2);
    end
    R12 = B; A2 = A(:,(kf+1):n);
end