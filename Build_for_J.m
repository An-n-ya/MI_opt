clc
clear
P = 2;
L = 5;
M = 3;
Q = 4;

for p = 1:P
    c(:,p) = rand(L,1);
    r(:,p) = rand(M,1);
    X_bar(:,1+(p-1)*M:p*M) = toeplitz(exp(2*pi*1i*c(:,p)),exp(2*pi*1i*r(:,p)));
end
X = kron(eye(Q),X_bar); 

m = exp(2*pi*1i*rand(M*L*P*Q^2,1));
x = vec(X');

J = zeros(M*L*Q*Q*P,(L+M-1)*P);
uni_x = zeros((L+M-1)*P,1);
s = 0;
for row = 1:(L+M-1)*P
    q = ceil(row/(L+M-1));
    i = mod(row,L+M-1);
    
    if i == 0
        i = L+M-1;
    end
    cur_x = X_bar(max(1,i-M+1),(q-1)*M+max(1,M-i+1));
    
    J(:,row) = vec(X' == cur_x');
    
%     if all(cur_x == x(col))
%         pause
% %         x(col)=0;
%     end
    uni_x(row,:) = cur_x;


end
 m'*x-m'*J*conj(uni_x);
 J'*J;
% uni_x'*J'*J*uni_x-x'*x
