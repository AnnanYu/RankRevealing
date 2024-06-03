function [Q,A] = givensr(A)
Q = eye(size(A,1));
for j=1:size(A,2)
  for i=j+1:size(A,1)
    rot = givens_rot(A(j,j),A(i,j));
    A([j,i],j:end) = rot' * A([j,i],j:end);
    Q(:,[j,i]) = Q(:,[j,i]) * rot;
  end
  A(j+1:end,j) = 0;
end

%only positive values on diagonal
D = diag(sign(diag(A)));
Q = Q*D; 
A = D*A;



