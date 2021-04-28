function [X,MI,MI_max] = Proj_MI(para)
    %% Initialization
    m = para.m;
    L = para.L;
    Q = para.Q;
    P = para.P;
    sigma = para.sigma;
    Sigma_g = para.Sigma_g;
    X0 = para.X0;
    P_max = para.P_max;
    
    %% Alternating Projection for PC
    t = 1;
    epsilon = 1e-3;
    delta = 1;
    X = X0;
    [U,D] = eig(Sigma_g);
    eta = water_filling(sigma^2./diag(D),P_max); 
    Omega_LHS = eta * eye(P*Q*m) - sigma^2 ./ D;
    Omega_LHS(Omega_LHS < 0) = 0.0001;
    Omega = sqrtm(Omega_LHS) * U';
    MI(t) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
    
    MI_t = log(sigma^(-2) * diag(D) * eta);
    MI_t(MI_t < 0) = 0;
    MI_max = sum(MI_t);
    
    clear X_bar;
    while abs(delta) > epsilon


        % Update Psi
        [V,~,W] = svd(Omega * X');
        if Q*P*m - Q*L > 0
            I = [eye(Q*L) zeros(Q*L, Q*P*m - Q*L)];
        else
            I = [eye(Q*P*m); zeros(Q*L - Q*P*m, Q*P*m)];
        end
        Psi = W * I * V';

        % Updaye X
        H = Psi * Omega;
        for k = 1:L
            for l = 1:P*m
                H_hat = H(k:L:L*Q,l:P*m:P*Q*m);
                X_bar(k,l) = trace(H_hat)/Q;
            end
        end
        X = kron(eye(Q),X_bar);


        % calculate MI
        MI(t+1) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
        delta = MI(t+1) - MI(t);
        t = t+1;
    end
end