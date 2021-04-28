clc
clear
P = 3;
L = 5;
M = 4;
Q = 7;

for p = 1:P
    c(:,p) = rand(L,1);
    r(:,p) = rand(M,1);
    X_bar(1+(p-1)*L:p*L,:) = toeplitz(exp(2*pi*1i*c(:,p)),exp(2*pi*1i*r(:,p)));
end
X = kron(eye(Q),X_bar); 

m = exp(2*pi*1i*rand(M*L*P*Q^2,1));
x = vec(X);

s = 0;
for row = 1:(L+M-1)*P
    q = ceil(row/(L+M-1));
    i = mod(row,L+M-1);
    
    if i == 0
        i = L+M-1;
    end
    cur_x = X_bar((q-1)*L+max(1,i-M+1),max(1,M-i+1));
    if i <= min(M,L)
        col = vec(M*L*P*Q^2-(L*(P-(q-1))-i) - (L*P*Q*(0:i-1)'+(0:i-1)')- M*L*P*Q*(0:Q-1)-L*P*max(0,Q-1)*(0:Q-1)/(Q-1));
    elseif L+M - i < min(M,L)
        col = vec(M*L*P*Q^2- L*P*Q*(i - L)-L*(p-q) - (L*P*Q*(0:L+M-i-1)'+(0:L+M-i-1)')- M*L*P*Q*(0:Q-1)-L*P*max(0,Q-1)*(0:Q-1)/(Q-1));
    elseif M >= L
        col = vec(M*L*P*Q^2- L*P*Q*(i - min(M,L))-L*(p-q) - (L*P*Q*(0:min(M,L)-1)'+(0:min(M,L)-1)')- M*L*P*Q*(0:Q-1)-L*P*max(0,Q-1)*(0:Q-1)/(Q-1));
    else
        col = vec(M*L*P*Q^2-(L*(P-(q-1))-i) - (L*P*Q*(0:min(M,L)-1)'+(0:min(M,L)-1)')- M*L*P*Q*(0:Q-1)-L*P*max(0,Q-1)*(0:Q-1)/(Q-1));
    end
%     if all(cur_x == x(col))
%         pause
% %         x(col)=0;
%     end
    uni_m(row) = sum(m(col));
    s = s+sum(conj(m(col)))*cur_x;
    uni_x(row) = cur_x;
    [repeat,~]=size(col);
    d(row) = repeat;
end
m'*x-s
uni_x*diag(d)*uni_x'-x'*x
