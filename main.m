%% Initialization
clc
clear
para.m = 5; %order of the FIR system
para.L = 50; %duration of the observed signals
para.Q = 2; %number of receive antenna
para.P = 2; %number of trasmit antenna
para.constant = 1; %constant modulus constraint
para.sigma = sqrt(2);
para.mode = 0; %normal mode
para.size1 = para.Q*para.P*para.m; %row number of calligraphic X
para.size2 = para.Q * para.L; %column number of calligraphic X
para.size3 = para.L+para.m-1; %number of distinct entries in each Xp
D = diag(rand(para.size1,1));
U = orth(rand(para.size1));
para.Sigma_g = U' * D * U;

build_sigma()

% construt the Toeplitz and Kronecker structure
for p = 1:para.P
    c(:,p) = rand(para.L,1);
    r(:,p) = rand(para.m,1);
    X_bar(:,1+(p-1)*para.m:p*para.m) = toeplitz(exp(2*pi*1i*c(:,p)),exp(2*pi*1i*r(:,p)));
end
para.X0 = kron(eye(para.Q),X_bar); 

J = zeros(para.size1 * para.size2,para.size3*para.P);

for row = 1:para.size3*para.P
    q = ceil(row/para.size3);
    i = mod(row,para.size3);
    if i == 0
        i = para.size3;
    end
    cur_x = X_bar(max(1,i-para.m+1),(q-1)*para.m+max(1,para.m-i+1));
    J(:,row) = vec(para.X0' == cur_x');
end
para.J = J;
para.D = J'*J;
% para.P_max = trace(para.X0'*para.X0) + 100; % power constraint
para.P_max = para.Q * para.L * 20;
%% Applying Algorithms

[X_acc_CMC,MI_acc_CMC] = MM_SQUAREM(para,'CMC');
[X_MM_CMC,MI_MM_CMC,~] = MM_MI(para,'CMC');


[X_acc_PC,MI_acc_PC] = MM_SQUAREM(para,'PC');
[X_MM,MI_MM,~] = MM_MI(para,'PC');
[X_proj,MI_Proj,MI_max] = Proj_MI(para);


%% Visualization

MSize = 6;
LWidth = 1.5;
figure(1)
plot(MI_acc_PC,'b-d','MarkerSize',MSize,'LineWidth',LWidth)
hold on
plot(MI_MM,'b--','MarkerSize',MSize,'LineWidth',LWidth)
hold on 
plot(MI_Proj,'r-x','MarkerSize',MSize,'LineWidth',LWidth)
legend('MM-SQUAREM','MM','Alt-Proj');
xlabel('Iterations');
ylabel('MI in nats');
print('PC20','-depsc','-painters')

% figure(2)
% plot(MI_acc_CMC,'b-d','MarkerSize',MSize,'LineWidth',LWidth)
% hold on 
% plot(MI_MM_CMC,'b--','MarkerSize',MSize,'LineWidth',LWidth)
% legend('MM-SQUAREM','MM');
% xlabel('Iterations');
% ylabel('MI in nats');
% print('CMC','-depsc','-painters')