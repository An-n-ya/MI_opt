N=8;
D = diag(rand(N,1));
U = orth(rand(N,N));

B = U' * D * U;
D = diag(rand(N,1));
U = orth(rand(N,N));
Bt = U' * D * U;


A = rand(N,N/2);
At = rand(N,N/2);

% LHS = log(det(eye(N/2) + A'*B*A));
% RHS = log(det(eye(N/2) + At'*B*At)) - trace(At'*B*At) + trace((A' * inv(inv(B) + A*A') *A)*(eye(N/2) + At'*B*At));
% LHS - RHS

LHS = trace(A'*inv(B)*A);;
RHS = 2 * real(trace(At'*inv(Bt)*A)) - trace(inv(Bt)*At*At'*inv(Bt)*B);
LHS - RHS