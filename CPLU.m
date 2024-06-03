function [Pi1,Pi2] = CPLU(A,kf)
    X = A; [m,n] = size(A);
    Pi1 = eye(m); Pi2 = eye(n);
    for k = 1:kf
        [maxc, rowindices] = max(abs(A(k:n, k:n)));
        [maxm, colindex] = max(maxc);
        row = rowindices(colindex) + k - 1; col = colindex + k - 1;
        A( [k, row], : ) = A( [row, k], : );
        A( :, [k, col] ) = A( :, [col, k] );
        Pi1( [k, row], : ) = Pi1( [row, k], : );
        Pi2( :, [k, col] ) = Pi2( :, [col, k] );
        A((k+1):m,k) = A((k+1):m,k) / A(k,k);
        i = k+1:n;
        A(i,i) = A(i,i) - A(i,k) * A(k,i);
    end
end