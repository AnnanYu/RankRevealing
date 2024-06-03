function [rat_lwb,rat_upb,i,j,s,t] = gammaLU(varargin)
% gammaLU    Compute the local maximum volume ratio in an LU factorization
% 
% [rat_lwb,rat_upb,i,j,s,t] = gammaLU(A,k) returns the value
%   of rat_lwb = max(vol(tilde(A_11)) / vol(A_11)), where A_11 = A(1:k,1:k)
%   and tilde(A_11) ranges over all k-by-k matrices
%   obtained from A_11 by a single row swap or a single column swap, or both.
%   In addition, the bound can be achieved
%   by swapping row i with row s, and column j with column t.
% 
% [rat_lwb,rat_upb,i,j,s,t] = gammaLU(A,k,tol) returns a lower bound
%   and an upper bound on max(vol(tilde(A_11)) / vol(A_11)), 
%   i.e., rat_lwb <= max(vol(tilde(A_11)) / vol(A_11)) <= rat_upb.
%   In addition, the lower bound can be achieved
%   by swapping row i with row s, and column j with column t.
%   The ratio between the upper bound and the lower bound 
%   is at most tol, i.e., rat_upb/rat_lwb < tol.

    A = varargin{1};
    k = varargin{2};
    if nargin == 2
        tol = 1;
    elseif nargin == 3
        tol = varargin{3};
    else
        error('maxvolratio accepts two or three parameters.\n');
    end

    [m,n] = size(A);
    if k == m
        A = [A;zeros(1,n)];
        [m,n] = size(A);
    end
    if k == n
        A = [A,zeros(m,1)];
        [m,n] = size(A);
    end

    % Compute auxiliary matrices
    A11inv = inv(A(1:k,1:k));
    U = A(1:k,1:k) \ A(1:k,(k+1):n);
    VT = A((k+1):m,1:k) / A(1:k,1:k);
    E = A((k+1):m,(k+1):n) - A((k+1):m,1:k) * U;

    % Sort auxiliary matrices
    A11inv_flat = flat_mat(A11inv,0,0); E_flat = flat_mat(E,k,k);
    VT_flat = flat_mat(VT,k,0); U_flat = flat_mat(U,0,k);

    % Compute the current upper and lower bounds
    Vmax = abs(VT_flat(1,1)); Umax = abs(U_flat(1,1));
    Amax = abs(A11inv_flat(1,1)); Emax = abs(E_flat(1,1));
    [rat_lwb,i,j,s,t] = update_lwb(A11inv_flat,E_flat,VT_flat,U_flat,k,...
        A11inv,E,VT,U,-1,[0,0,0,0]);
    rat_upb = max([Vmax, Umax, Vmax*Umax + Amax*Emax, rat_lwb]);

    % Test rank-two updates until tolerance reached
    iter = 0; MAXIT = 1e4; warning_flag = true;
    while rat_upb / rat_lwb > tol
        % Test all cases based on swapping the ssth row in
        ss = VT_flat(1,2); inds_col = find(VT_flat(:,2) == ss);
        inds_col_E = E_flat(:,2) == ss;
        for icol = 1:length(inds_col)
            ii = VT_flat(inds_col(icol),3); v = VT_flat(inds_col(icol),1);
            [rat_lwb,i,j,s,t] = sweep_lower(rat_lwb,v,U_flat(:,1),...
                A11inv(U_flat(:,2)',ii),transpose(E(ss-k,U_flat(:,3)'-k)),...
                ii,ss,[i,j,s,t],U_flat); iter = iter + 1;
        end

        % Test all cases based on swapping the iith row out
        ii = VT_flat(1,3); inds_row = find(VT_flat(:,3) == ii);
        inds_row_A = A11inv_flat(:,3) == ii;
        for irow = 1:length(inds_row)
            ss = VT_flat(inds_row(irow),2); v = VT_flat(inds_row(irow),1);
            [rat_lwb,i,j,s,t] = sweep_lower(rat_lwb,v,U_flat(:,1),...
                A11inv(U_flat(:,2)',ii),transpose(E(ss-k,U_flat(:,3)'-k)),...
                ii,ss,[i,j,s,t],U_flat); iter = iter + 1;
        end

        % Remove entries that have been considered
        VT_flat(union(inds_col,inds_row),:) = [];
        E_flat(inds_col_E,:) = [];
        A11inv_flat(inds_row_A,:) = [];

        if iter > MAXIT && warning_flag
            warning(strcat('maxvolratio has not reached tolerance.\n',...
            'Current rat_upb/rat_lwb: %.3e'), rat_upb / rat_lwb);
            warning_flag = false;
        end

        if ~isempty(VT_flat) && ~isempty(U_flat)
            % Recompute the lower and upper bounds
            % Uncommenting this line does not seem to change a lot
            [rat_lwb,i,j,s,t] = update_lwb(A11inv_flat,E_flat,VT_flat,U_flat,k,...
                A11inv,E,VT,U,rat_lwb,[i,j,s,t]);
            rat_upb = max([Vmax, Umax, rat_lwb, ...
                abs(VT_flat(1,1)) * abs(U_flat(1,1)) ...
                + abs(A11inv_flat(1,1)) * abs(E_flat(1,1))]);
        else
            % All swaps considered.
            rat_upb = rat_lwb; break
        end
    end
end

function [curr_lower,i,j,s,t] = update_lwb(A11inv_flat,E_flat,VT_flat,U_flat,k,...
    A11inv,E,VT,U,curr_lower,inds)
    i = inds(1); j = inds(2); s = inds(3); t = inds(4);
    Vmax = abs(VT_flat(1,1)); Umax = abs(U_flat(1,1));
    Amax = abs(A11inv_flat(1,1)); Emax = abs(E_flat(1,1));
    [curr_lower,maxind] = max([Vmax, Umax, Vmax*Umax - Amax*Emax, ...
        -Vmax*Umax + Amax*Emax, 1, curr_lower]);
    switch maxind
        case 1
            i = VT_flat(1,3); j = 1; s = VT_flat(1,2); t = 1;
        case 2
            i = 1; j = U_flat(1,2); s = 1; t = U_flat(1,3);
        case 3
            i = VT_flat(1,3); j = U_flat(1,2);
            s = VT_flat(1,2); t = U_flat(1,3);
            curr_lower = abs(VT(s-k,i)*U(j,t-k) + A11inv(j,i)*E(s-k,t-k));
        case 4
            i = A11inv_flat(1,3); j = A11inv_flat(1,2);
            s = E_flat(1,2); t = E_flat(1,3);
            curr_lower = abs(VT(s-k,i)*U(j,t-k) + A11inv(j,i)*E(s-k,t-k));
        case 5
            i = 1; j = 1; s = 1; t = 1;
    end
end

function A_flat = flat_mat(A,augx,augy)
    [m,n] = size(A);
    [rowind,colind] = meshgrid(1:m,1:n);
    A_flat = [reshape(A,[],1),reshape(rowind',[],1),reshape(colind',[],1)];
    [~,sortind] = sort(abs(A_flat(:,1)),'descend');
    A_flat = A_flat(sortind,:);
    A_flat(:,2) = A_flat(:,2) + augx; A_flat(:,3) = A_flat(:,3) + augy;
end

function [curr_lower,i,j,s,t] = sweep_lower(curr_lower,v,u,a,e,ii,ss,inds,U_flat)
    i = inds(1); j = inds(2); s = inds(3); t = inds(4);
    Rats = v * u + a .* e; [rat,maxind] = max(abs(Rats));
    if rat > curr_lower
        curr_lower = rat;
        i = ii; j = U_flat(maxind,2); s = ss; t = U_flat(maxind,3);
    end
end