global Q; Q = para.Q;
global P; P = para.P;
global m; m = para.m;
global LR; LR = para.LR;
global LT; LT = para.LT;
size1 = P*Q*m;
% D = diag(rand(size1,1));
% U = orth(rand(size1));
global Sigma_g; Sigma_g = U' * D * U;
% g = mvnrnd(zeros(size1),Sigma_g);

OFFSET = 1;
SIGMA = zeros(P * Q * LT * LR);
for x = 1 : LR * Q
    for y = 1 : LR * Q
        SIGMA((x-1)*LT*P+1:x*LT*P, (y-1)*LT*P+1:y*LT*P) = covariance(x,y,OFFSET);
    end
end

clear D and m and P and Q and size1 and U and x and y and OFFSET and Sigma_g

function g_entry = find_g(i,j,t1,p,q,t2)
    global P;
    global m;
    global Sigma_g;
    if t1 > m || t1 <= 0 || t2 > m || t2 <= 0
        g_entry = 0;
    else
        g_x = P * m * j + m * i + t1;
        g_y = P * m * q + m * p + t2;
        g_entry = Sigma_g(g_x,g_y);
    end
end

function c = covariance(u,k,bias)
    global LR;
    global LT;
    global P;
    
    c = zeros(P*LT);
    q = floor(u / LR);
    j = floor(k / LR);
    
    bias1 = mod(u,LR);
    if bias1 == 0
        q = q - 1;
        bias1 = LR;
    end
    t2 = bias + bias1;
    
    
    bias2 = mod(k,LR);
    if bias2 == 0
        j = j - 1;
        bias2 = LR;
    end
    t1 = bias + bias2;
    
    for iter1 = 1:P*LT
        p = floor(iter1 / LT);
        bias11 = mod(iter1,LT)-1;
        if bias11 < 0
            p = p - 1;
            bias11 = LT-1;
        end
        for iter2 =1:P*LT
            i = floor(iter2 / LT);
            bias22 = mod(iter2,LT)-1;
            if bias22 < 0
                i = i - 1;
                bias22 = LT-1;
            end
            c(iter1,iter2) = find_g(i,j,t1-bias22,p,q,t2-bias11);
        end
    end
end