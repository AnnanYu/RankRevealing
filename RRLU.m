function [L,U,Pi1,Pi2,Ak,swaps,t1,t2] = RRLU(A,kf,rho)
% RRLU    Compute a strong rank-revealing LU factorization
% 
% [L,U,Pi1,Pi2,Ak,swaps] = RRLU(A,kf,rho) returns a rank-kf partial
% LU factorization of a matrix A so that the pivots Pi1(1:kf,:)*A*Pi2(:,1:kf)
% satisfies a near local maximum volume property with a ratio
% bounded by rho. The outputs satisfy that
% L*U = Pi1(1:kf,:)*A*Pi2(:,1:kf) and that Ak is the
% rank-kf approximation of Pi1*A*Pi2 constructed by Gaussian
% elimination. The output swaps contains the number of swaps being
% performed to obtain a near local maximum subset. It is recommended that
% rho >= 3.

    tic
    rho = rho / 3 * 2;
    [m,n] = size(A); swaps = 0; k = kf;
    if k > min([m,n])
        throw(MException('InputError:rankOverflow',...
            'k cannot be larger than m or n'));
    end
    if rho < 2
        throw(MException('InputError:rhoTooSmall',...
            'rho should not be smaller than 3'));
    end
    [Pi1,Pi2] = CPLU(A,kf); % GE with complete pivoting
    A = Pi1 * A * Pi2;
    t1 = toc;

    if k == min([m,n])
        Ak = Pi1 * A * Pi2; [L,U] = lu(Ak); t2 = toc; return
    end

    % Compute the auxiliary matrices needed in the rank update formulas
    A11 = A(1:k,1:k); A11inv = inv(A11); U = A((k+1):end,1:k) / A11;
    VT = A11 \ A(1:k,(k+1):end);
    E22 = A((k+1):end,(k+1):end) - A((k+1):end,1:k) * VT;
    subcount = 0; MAXIT = 100;
    while subcount <= MAXIT
        subcount_inner = 0;
        % Do all rank-1 updates
        while (~isempty(VT) && max(abs(VT),[],'all') > sqrt(rho/2)) || ...
                (~isempty(U) && max(abs(U),[],'all') > sqrt(rho/2))
            subcount_inner = subcount_inner + 1;
            % If the while-loop runs for over 20 iterations, then a
            % numerical error has likely occured. In that case, recompute
            % the auxiliary matrices
            if subcount_inner > 20
                [U,VT,A11inv,E22] = recompute(A,k); subcount_inner = 0;
            end
            swaps = swaps + 1;
            if max(abs(VT),[],'all') > max(abs(U),[],'all')
                % Swap a column
                [~,I] = max(abs(VT),[],'all'); [s,t] = ind2sub(size(VT),I);
                Pi2(:,[s,k+t]) = Pi2(:,[k+t,s]);

                % Here, update the four auxiliary matrices
                L = A(:,k+t) - A(:,s);
                L1 = L(1:k,:); L2 = L((k+1):m,:);
                R = zeros(1,n); R(s) = 1; R(k+t) = -1;
                R1 = R(:,1:k); R2 = R(:,(k+1):n);
                [U,VT,A11inv,E22] = update_aux(U,VT,A11inv,E22,L1,L2,R1,R2);
                A(:,[s,k+t]) = A(:,[k+t,s]);
            else
                % Swap a row
                [~,I] = max(abs(U),[],'all'); [j,i] = ind2sub(size(U),I);
                Pi1([i,k+j],:) = Pi1([k+j,i],:);

                % Here, update the four auxiliary matrices
                L = zeros(m,1); L(i) = 1; L(k+j) = -1;
                L1 = L(1:k,:); L2 = L((k+1):m,:);
                R = A(k+j,:) - A(i,:);
                R1 = R(:,1:k); R2 = R(:,(k+1):n);
                [U,VT,A11inv,E22] = update_aux(U,VT,A11inv,E22,L1,L2,R1,R2);
                A([i,k+j],:) = A([k+j,i],:);
            end
        end

        % Check for a rank-2 update
        [M1,I1] = max(abs(A11inv),[],'all');
        [M2,I2] = max(abs(E22),[],'all');
        if ~isempty(M2) && M1*M2 > rho
            subcount = subcount + 1;
            if mod(subcount,20) == 0
                [U,VT,A11inv,E22] = recompute(A,k);
            end
            swaps = swaps + 1;
            [s,i] = ind2sub(size(A11inv),I1);
            [j,t] = ind2sub(size(E22),I2);
            Pi1([i,k+j],:) = Pi1([k+j,i],:);
            Pi2(:,[s,k+t]) = Pi2(:,[k+t,s]);

            % Here, update the four auxiliary matrices;
            L = [A(:,k+t) - A(:,s), zeros(m,1)]; L(i,2) = 1; L(k+j,2) = -1;
            L1 = L(1:k,:); L2 = L((k+1):m,:);
            R = [A(k+j,:) - A(i,:); zeros(1,n)]; R(2,s) = 1; R(2,k+t) = -1;
            R1 = R(:,1:k); R2 = R(:,(k+1):n);
            [U,VT,A11inv,E22] = update_aux(U,VT,A11inv,E22,L1,L2,R1,R2);
            A([i,k+j],:) = A([k+j,i],:);
            A(:,[s,k+t]) = A(:,[k+t,s]);
        else
            break
        end
    end
    if subcount > MAXIT
        warning(['The requested rho has not been reached after 100 swaps.' ...
            'The algorithm terminates without a guarantee on rho']);
    end
    Ak = [eye(kf);U] * A(1:kf,1:kf) * [eye(kf),VT];
    [L,U] = lu(A(1:kf,1:kf));
    t2 = toc;
end