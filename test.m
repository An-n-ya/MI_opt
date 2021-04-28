%% Initialization
clc
clear
m = 4;
L = 3;
Q = 2;
P = 5;
constant = 1;
sigma = 1;
% n = random('Normal',0,sqrt(sigma/2),Q*P*M,1)+1j*random('Normal',0,sqrt(sigma/2),Q*P*M,1);

D = diag(rand(Q*P*m,1));
U = orth(rand(Q*P*m));
Sigma_g = U' * D * U;

for p = 1:P
    c(:,p) = rand(L,1);
    r(:,p) = rand(m,1);
    X_bar(:,1+(p-1)*m:p*m) = toeplitz(exp(2*pi*1i*c(:,p)),exp(2*pi*1i*r(:,p)));
end
X0 = kron(eye(Q),X_bar); 

J = zeros(m*L*Q*Q*P,(L+m-1)*P);
for row = 1:(L+m-1)*P
    q = ceil(row/(L+m-1));
    i = mod(row,L+m-1);  
    if i == 0
        i = L+m-1;
    end
    cur_x = X_bar(max(1,i-m+1),(q-1)*m+max(1,m-i+1));   
    J(:,row) = vec(X0' == cur_x');
end

%% test
t = 1;
epsilon = 1e-7;
delta = 1;
X = X0;
MI(t) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));

    X_old = X;
    M_LHS = eye(Q*L) + X * sigma^(-2) * Sigma_g * X';
    M_RHS = X / (sigma^2 * inv(Sigma_g) + X' * X);
    M = M_LHS * M_RHS;
    N = M_RHS' * M_LHS * M_RHS;
    PP = kron(eye(Q*L),N);
    lambda_P = max(real(eig(N)));
    z = J' * (vec(M') - (PP - lambda_P * eye(Q*Q*L*m*P))* vec(X'));
    x_hat = conj(constant * exp(1j * angle(z)));
    X = (reshape(J * x_hat,[Q*m*P, Q*L]))';
    MI(t+1) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
    delta = MI(t+1) - MI(t);
    t = t+1;

    Xt = X_old;
    RHS1 = log(det(eye(6) + Xt*Sigma_g*Xt')) - trace(Xt*Sigma_g*Xt') + trace((X * inv(inv(Sigma_g) + X'*X) *X')*(eye(6) + Xt*Sigma_g*Xt'));
    LHS1 = log(det(eye(6)+X*Sigma_g*X'));
    res1 = LHS1 - RHS1
    
    LHS2 = trace((X * inv(inv(Sigma_g) + X'*X) *X')*(eye(6) + Xt*Sigma_g*Xt'));
    LHS22 = trace(sqrtm(eye(6) + Xt*Sigma_g*Xt')*X*inv(inv(Sigma_g) + X'*X)*X'*sqrtm(eye(6) + Xt*Sigma_g*Xt'));
    B = inv(Sigma_g) + X'*X;
    A = X'*sqrtm(eye(6) + Xt*Sigma_g*Xt');
    At = Xt'*sqrtm(eye(6) + Xt*Sigma_g*Xt');
    Bt = inv(Sigma_g) + Xt'*Xt;
    res = trace(A'*inv(B)*A)-2 * real(trace(At'*inv(Bt)*A)) + trace(inv(Bt)*At*At'*inv(Bt)*B);
    
    RHS2 = 2 * real(trace(M * X')) - trace(N * (inv(Sigma_g) + X'*X));
    res2 = LHS2 - RHS2
    
    
    
    x = vec(X');
    mm = vec(M');
    xt = vec(X_old');
    LHS3 = 2*real(mm'*x)-x'*PP*x;
    RHS3 = 2 * real(mm'*x) - lambda_P * x'*x - 2 * real(x' * (PP - lambda_P * eye(240)) * xt) + xt'*(lambda_P * eye(240) - PP)*xt;
    LHS3 - RHS3
    
    RHS33 = 2 * real((mm-(PP - lambda_P * eye(240)) * xt)'*x) - lambda_P * x'*x + xt'*(lambda_P * eye(240) - PP)*xt;
    
%     err_L = trace(A'*inv(B)*A) - LHS2
%     err_R = 2 * real(trace(At'*inv(Bt)*A)) - trace(inv(Bt)*At*At'*inv(Bt)*B) - RHS2
    