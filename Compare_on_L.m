%% Initialization
clc
clear
para.m = 5; %order of the FIR system %when m gets larger, problem happens too.
para.L = 4; %duration of the observed signals
para.Q = 2; %number of receive antenna
para.P = 2; %number of trasmit antenna %when P or Q is too large in terms of L, Proj algorithm will have problem.
para.constant = 1; %constant modulus constraint
para.sigma = sqrt(2);
para.mode = 0; %normal mode
para.size1 = para.Q*para.P*para.m; %row number of calligraphic X

D = diag(rand(para.size1,1));
U = orth(rand(para.size1));
para.Sigma_g = U' * D * U;

%% Comparison
ind = 1;
for i = 2:2:16
    %sub-initial
    para.L = i;
    para.size2 = para.Q * para.L; %column number of calligraphic X
    para.size3 = para.L+para.m-1; %number of distinct entries in each Xp
    clear c r X_bar para.X0 J para.J para.D
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
    [X_acc_PC,MI_acc_PC] = MM_SQUAREM(para,'PC');
    [X_MM,MI_MM,~] = MM_MI(para,'PC');
    [X_proj,MI_Proj,MI_max] = Proj_MI(para);
    mi_acc_PC(ind) = MI_acc_PC(end);
    mi_MM_PC(ind) = MI_MM(end);
    mi_Proj(ind) = MI_Proj(end);
    mi_max(ind) = MI_max;
    ind = ind + 1;
end


%% Visualization
MSize = 6;
LWidth = 1.5;
figure(1)
plot(2:2:16,mi_acc_PC,'b-d','MarkerSize',MSize,'LineWidth',LWidth)
hold on 
plot(2:2:16,mi_MM_PC,'b-o','MarkerSize',MSize,'LineWidth',LWidth)
hold on 
plot(2:2:16,mi_Proj,'r-x','MarkerSize',MSize,'LineWidth',LWidth)
hold on 
plot(2:2:16,mi_max,'k-*','MarkerSize',MSize,'LineWidth',LWidth)
legend('MM-SQUAREM','MM','Alt-Proj','Analytical-maximum');
xlabel('Number of L');
ylabel('MI in nats');
% print('Comparison_on_L','-depsc','-painters')