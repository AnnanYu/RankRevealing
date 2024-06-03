function [Q,R,Pi,Ak,swaps,t1,t2] = RRQR(A,kf,rho)
% RRQR    Compute a strong rank-revealing QR factorization
% 
% [Q,R,Pi,Ak,swaps] = RRQR(A,kf,rho) returns a rank-kf partial
% QR factorization of a matrix A so that the pivots A*Pi(:,1:kf)
% satisfies a near local maximum volume property with a ratio
% bounded by rho. The outputs satisfy that Q*R = A*Pi(:,1:kf) and that
% Ak is the rank-kf approximation of A*Pi constructed by partial QR.
% The output swaps contains the number of swaps being
% performed to obtain a near local maximum subset. It is recommended that
% rho >= 2.

    tic
    rho = rho / sqrt(2);
    [m,n] = size(A); k = kf; swaps = 0;
    if k > min([m,n])
        throw(MException('InputError:rankOverflow',...
            'k cannot be larger than m or n'));
    end
    if rho < 1
        throw(MException('InputError:rhoTooSmall',...
            'rho should not be smaller than sqrt(2)'));
    end
    [Q,R,Pi,gamma,R12,A2] = CPQR(A,kf); % QR with column pivoting
    t1 = toc;

    if k == min([m,n])
        Ak = A * Pi; t2 = toc; return
    end
    
    % Compute the auxiliary matrices
    gamma = gamma((k+1):end); R11inv = inv(R);
    R11invR12 = R \ R12; omega = vecnorm(transpose(R11inv)).^-1;

    % Check for local maximum volume
    [rhohat_eval,swt] = max([max(max(abs(R11invR12))),max(gamma)/min(omega)]);

    while rhohat_eval > rho
        if swt == 1
            [M,I] = max(abs(R11invR12));
            [~,ind_c] = max(M); ind_r = I(ind_c);
        else
            [~,ind_c] = max(gamma); [~,ind_r] = min(omega); 
        end

        %prod(svd(A*Pi(:,1:5)))
        swaps = swaps + 1;

        % Reduce to the case where we swap column k and k+1
        R12(:,[1,ind_c]) = R12(:,[ind_c,1]);
        Pi(:,[k+1,k+ind_c]) = Pi(:,[k+ind_c,k+1]);
        R11invR12(:,[1,ind_c]) = R11invR12(:,[ind_c,1]);
        A2(:,[1,ind_c]) = A2(:,[ind_c,1]);
        gamma([1,ind_c]) = gamma([ind_c,1]);

        if ind_r ~= k
            Pi(:,[ind_r,k]) = Pi(:,[k,ind_r]);
            R(:,[ind_r,k]) = R(:,[k,ind_r]);
            R11invR12([ind_r,k],:) = R11invR12([k,ind_r],:);
            omega([ind_r,k]) = omega([k,ind_r]);

            % reformat with givens (yet givens here inefficient)
            [Q_giv,R1] = givensr(R(ind_r:k, ind_r:k));
            R(ind_r:k, ind_r:k) = R1;
            R12(ind_r:k, :) = Q_giv' * R12(ind_r:k, :);
            R11inv(:, ind_r:k) = R11inv(:, ind_r:k) * Q_giv;
            R11inv([ind_r,k],:) = R11inv([k,ind_r],:);
            Q(:,ind_r:k) = Q(:,ind_r:k) * Q_giv;
        end

        % Update the auxiliary matrices
        Pi(:,[k,k+1]) = Pi(:,[k+1,k]);
        A2 = A2 + Q(:,k) * R12(k,:); % Backtrack
        a_k = Q(:,k)*R(k,k);
        q = A2(:,1); nq = norm(q); q = q/nq; Q = [Q(:,1:(k-1)),q]; % new kth qr column
        Rnew = [R(1:(k-1),1:(k-1)),R12(1:(k-1),1);zeros(1,k-1),nq];
        b = q' * a_k; % Project a_k onto q
        Aonq = q' * A2(:,2:end); % Project A2 onto q
        R12new = [[R(1:(k-1),k);b],[R12(1:(k-1),2:end); Aonq]]; % Compute new R12
        A2_new = [a_k - q*b, A2(:,2:end) - q * Aonq]; % Compute new residual
        u = Rnew(:,k) - R(:,k); v = zeros(1,k); v(k) = 1; % Find the rank-1 update of R11
        x = zeros(k,1); x(k) = 1; y = R12new(k,:) - R12(k,:); % Find the rank-1 update of R12
        Rinvu = R \ u; vRinv = v / R; vRinvu = v * Rinvu;
        Rinvnew = R11inv - Rinvu * vRinv / (1+vRinvu); % Use Woodbury to recompute R11^-1
        RinvRnew = R11invR12 + (R \ x) * y - Rinvu * (vRinv * R12new) / (1+vRinvu); % Use Woodbury to recompute R11^-1 R12
        RinvRnew(:,1) = Rnew \ R12new(:,1);

        R = Rnew; R12 = R12new; A2 = A2_new; gamma = vecnorm(A2);
        R11inv = Rinvnew; R11invR12 = RinvRnew; omega = vecnorm(transpose(R11inv)).^-1;

        [rhohat_eval,swt] = max([max(max(abs(R11invR12))),max(gamma)/min(omega)]);
    end
    Ak = Q * [R R12];
    t2 = toc;
end