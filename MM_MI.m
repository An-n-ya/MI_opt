function [X,MI,x] = MM_MI(para, mode)
    %% Initialization
    m = para.m;
    LT = para.LT;
    LR = para.LR;
    Q = para.Q;
    P = para.P;
    sigma = para.sigma;
    Sigma_g = para.Sigma_g;
    X0 = para.X0;
    x0 = para.x0;
    SIGMA = para.SIGMA;
    inv_SIGMA = para.inv_SIGMA;
    
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
        x = x0;
        MI(t) = real(log(det(sigma^(-2)*X'*SIGMA*X+eye(LR*Q))));
        while abs(delta) > epsilon
            V_LHS = eye(Q*LR) + X' * sigma^(-2) * SIGMA * X;
            V_RHS = X' / (sigma^2 * inv_SIGMA + X * X');
            V = V_LHS * V_RHS;
            U = V' * V_RHS;
            %vectorize V and extract W
            v = zeros(P*LT,1);
            W = zeros(P*LT);
            for i = 1:Q*LR
                ind_begin = (i-1)*LT*P+1;
                ind_end = i*LT*P;
                v = v + V(i,ind_begin:ind_end)';
                W = W + U(ind_begin:ind_end,ind_begin:ind_end);
            end
            n = LT*P;
            W(eye(n)==1) = real(W(eye(n)==1));
            l = chol(W);
            
            cvx_begin quiet
                variable x(n) complex
                maximize( 2*real(v'*x) - x'*(l'*l)*x )
                subject to
                    x'*x <= P_max
            cvx_end
            X = kron(eye(Q*LR),x); 
            MI(t+1) = real(log(det(sigma^(-2)*X'*SIGMA*X+eye(LR*Q))));
            delta = MI(t+1) - MI(t);
            t = t+1;
            if para.mode == 1
                delta = 0;
            end
        end
    elseif strcmp(mode, 'CMC')
           %% MM Algorithm for CMC
            t = 1;
            epsilon = 1e-4;
            delta = 1;
            X = X0;
            x = x0;
            MI(t) = real(log(det(sigma^(-2)*X'*SIGMA*X+eye(LR*Q))));
            while abs(delta) > epsilon
                V_LHS = eye(Q*LR) + X' * sigma^(-2) * SIGMA * X;
                V_RHS = X' / (sigma^2 * inv_SIGMA + X * X');
                V = V_LHS * V_RHS;
                U = V' * V_RHS;
                %vectorize V and extract W
                v = zeros(P*LT,1);
                W = zeros(P*LT);
                for i = 1:Q*LR
                    ind_begin = (i-1)*LT*P+1;
                    ind_end = i*LT*P;
                    v = v + V(i,ind_begin:ind_end)';
                    W = W + U(ind_begin:ind_end,ind_begin:ind_end);
                end
                lambda_W = max(real(eig(W)));
                
                %update x
                z = v - (W - lambda_W * eye(LT*P))* x;
                x = constant * exp(1j * angle(z));
                X = kron(eye(Q*LR),x); 
                MI(t+1) = real(log(det(sigma^(-2)*X'*SIGMA*X+eye(LR*Q))));
                delta = MI(t+1) - MI(t);
                t = t+1;
                if para.mode == 1
                    delta = 0;
                end
            end
    end
end

