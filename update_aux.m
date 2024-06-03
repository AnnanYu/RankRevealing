function [U,VT,A11inv,E22] = update_aux(U,VT,A11inv,E22,L1,L2,R1,R2)
    [~,d] = size(L1);
    RL11 = R1*A11inv*L1; L2RL11 = L2*RL11; RL11R2 = RL11*R2;
    A11invL1 = A11inv*L1; R1A11inv = R1*A11inv;
    UL1 = U*L1; R1VT = R1*VT;
    V = inv(eye(d) + RL11);
    A11inv_new = A11inv - A11invL1*V*R1A11inv;
    U_new = U + L2*R1A11inv - UL1*V*R1A11inv - L2RL11*V*R1A11inv;
    VT_new = VT + A11invL1*R2 - A11invL1*V*R1VT - A11invL1*V*RL11R2;
    E22_new = E22 + L2*R2 - L2*R1VT - UL1*R2 - L2RL11*R2 ...
        + L2RL11*V*R1VT + UL1*V*RL11R2 + UL1*V*R1VT + L2RL11*V*RL11R2;
    A11inv = A11inv_new; U = U_new; VT = VT_new; E22 = E22_new;
end

% Test script
% A11 = rand(5,5); A12 = rand(5,3); A21 = rand(3,5); A22 = rand(3,3);
% L1 = rand(5,2); L2 = rand(3,2); R1 = rand(2,5); R2 = rand(2,3);
% A11inv = inv(A11); U = A21 * A11inv; VT = A11inv * A12; E22 = A22 - U * A12;
% [U,VT,A11inv,E22] = update_aux(U,VT,A11inv,E22,L1,L2,R1,R2)
% A11 = A11 + L1 * R1; A12 = A12 + L1 * R2; A21 = A21 + L2 * R1; A22 = A22 + L2 * R2;
% A11inv = inv(A11)
% U = A21 * A11inv
% VT = A11inv * A12
% E22 = A22 - U * A12