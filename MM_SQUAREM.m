
function [X,MI] = MM_SQUAREM(para, mode)
    %% Initialization
    global m L Q P J sigma Sigma_g constant P_max D
    m = para.m;
    L = para.L;
    Q = para.Q;
    P = para.P;
    J = para.J;
    sigma = para.sigma;
    Sigma_g = para.Sigma_g;
    constant = para.constant;
    P_max = para.P_max;
    D = J'*J;
    X0 = para.X0;
    para.mode = 1; %single loop
    
    
    %% SQUAREM Acceleration
    t = 1;
    epsilon = 1e-3;
    delta = 1;
    X = X0;
    MI(t) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
    while abs(delta) > epsilon
        para.X0 = X;
        [~,~,x1] = MM_MI(para,mode);
        X1 = (reshape(J * conj(x1),[Q*m*P, Q*L]))';
        para.X0 = X1;
        [~,~,x2] = MM_MI(para,mode);
        x = extract(X);
        q = x1 - x;
        v = x2 - x1 - q;
        alpha =  -norm(q)/norm(v);
        x_new =  proj((x - 2 * alpha * q + alpha^2 * v),mode);
        cal_MI(x_new)
        while cal_MI(x_new) < cal_MI(x)
            cal_MI(x_new)
            alpha = (alpha - 1) / 2;
            x_new =  proj(x - 2 * alpha * q + alpha^2 * v,mode);
        end
        X = (reshape(J * x_new,[Q*m*P, Q*L]))';
        MI(t+1) = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
        delta = MI(t+1) - MI(t);
        t = t+1;
    end
    plot(MI)
end

function x = extract(X)
    global m L Q P
    [r,n] = size(X);
    X_bar = X(1:r/Q,1:n/Q);
    x = zeros((L+m-1)*P,1);
    for row = 1:(L+m-1)*P
        q = ceil(row/(L+m-1));
        i = mod(row,(L+m-1));  
        if i == 0
            i = (L+m-1);
        end
        cur_x = X_bar(max(1,i-m+1),(q-1)*m+max(1,m-i+1));   
        x(row) = cur_x;
    end
end

function x_proj = proj(xx,mode) 
    global m L Q P J sigma Sigma_g constant P_max D
    if strcmp(mode, 'PC')
        X = (reshape(J * conj(xx),[Q*m*P, Q*L]))';
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
        x_proj = x;
    elseif strcmp(mode, 'CMC')
        x_proj = conj(constant * exp(1j * angle(xx)));
    end
end

function mi = cal_MI(x)
    global m L Q P J sigma Sigma_g
    X = (reshape(J * conj(x),[Q*m*P, Q*L]))';
    mi = real(log(det(sigma^(-2)*Sigma_g*(X'*X)+eye(P*Q*m))));
end
