function [X,MI,x_hat] = MM_MI(para, mode)
    %% Initialization
    m = para.m;
    L = para.L;
    Q = para.Q;
    P = para.P;
    sigma = para.sigma;
    Sigma_g = para.Sigma_g;
    X0 = para.X0;
    
    J = para.J; 
    D = para.D;
    if strcmp(mode, 'PC')
        P_max = para.P_max;
    elseif strcmp(mode, 'CMC')
        constant = para.constant;
    end
    
    if strcmp(mode, 'PC')
        %% MM Algorithm for PC
        t = 1;
        epsilon = 1e-3;
        delta = 1;
        X = X0;
        MI(t) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
        while abs(delta) > epsilon
            M_LHS = eye(Q*L) + X * sigma^(-2) * Sigma_g * X';
            M_RHS = X / (sigma^2 * inv(Sigma_g) + X' * X);
            M = M_LHS * M_RHS;
            mm = vec(M');
            N = M_RHS' * M_LHS * M_RHS;
            PP = kron(eye(Q*L),N);
            [~,n] = size(J);
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
            if para.mode == 1
                delta = 0;
            end
        end
        x_hat = x;
    elseif strcmp(mode, 'CMC')
           %% MM Algorithm for CMC
            t = 1;
            epsilon = 1e-4;
            delta = 1;
            X = X0;
            MI(t) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
            while abs(delta) > epsilon
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
                if para.mode == 1
                    delta = 0;
                end
            end
    end
end

