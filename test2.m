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
D = J'*J;
P_max = trace(X0'*X0) + 100;

%% test
t = 1;
epsilon = 1e-4;
delta = 1;
X = X0;
MI(t) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));

    X_old = X;
    M_LHS = eye(Q*L) + X * sigma^(-2) * Sigma_g * X';
    M_RHS = X / (sigma^2 * inv(Sigma_g) + X' * X);
    M = M_LHS * M_RHS;
    mm = vec(M');
    N = M_RHS' * M_LHS * M_RHS;
    PP = kron(eye(Q*L),N);
    
%     Power_temp = conj(J'*PP*J);
%     x_t = inv(Power_temp) * conj(J' * mm);
%     lambda_max = 5;
%     lambda_min = 0;
%     if real(x_t'* D * x_t) <= P_max
%         x_hat = x_t;
%     else
% %         [V,Diag] = eig(Power_temp);
% 
%         P_cur = P_max + 50;
%         threshold = 10;
%         lambda = lambda_max;
%         while P_max - P_cur > threshold || P_max < P_cur
%             x_t = inv(Power_temp - lambda*D) * conj(J' * mm);
%             P_cur = x_t'*D*x_t;
%             if P_cur <= P_max
%                 lambda_max = lambda;
%                 lambda = (lambda_min + lambda_max) / 2;
%             elseif P_cur > P_max
%                 lambda_min = lambda;
%                 lambda = (lambda_min + lambda_max) / 2;
%             end
%         end
%         x_hat = x_t;
%     end 
    n = L*P*Q;
    U = J'*PP*J;
    U(eye(n)==1) = real(U(eye(n)==1));
    l = chol(U);
    cvx_begin quiet
        variable x(n) complex
        maximize( 2*real(mm'*J*conj(x)) - (conj(x))'*(l'*l)*conj(x) )
        subject to
            x'*D*x <= P_max
    cvx_end
    X = (reshape(J * conj(x),[Q*m*P, Q*L]))';
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
    
    
    
    2*real(trace(M * X')) - trace(X*N*X')
    2*real(trace(M * X0')) - trace(X0*N*X0')
    x = vec(X');
    xt = rand(30,1);
%     norm(2*real(conj(J'*mm)) - 2*real(Power_temp)*xt+ 2 *lambda * D * xt)
%     norm(2*real(conj(J'*mm)) - 2*real(Power_temp)*x_hat+ 2 *lambda * D * x_hat)
%     2 * real(mm' * J * conj(x_hat)) - conj(x_hat)'*J'*PP*J*conj(x_hat)
    f(x,J'*PP*J,J'*mm,lambda*D)
    df(x,J'*PP*J,J'*mm,lambda*D)
	function res = f(x,Z,m,D)
        res = -conj(x)'*Z*conj(x)+2*real(m'*conj(x))+conj(x)'*D*conj(x);
    end
    
    function res = df(x,Z,m,D)
        res = norm(-2*conj(Z)*x+2*conj(m)+2*conj(D)*x);
    end
    