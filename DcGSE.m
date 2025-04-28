function [U] = DcGSE(in_X,true_labs,alpha,lambda,delta)
%% size
max_iter = 30;
KMAXiter = 1000;
KREPlic = 10;
num_smp = size(in_X,2);
num_cluster = max(true_labs);
num_view = 3;
k = num_cluster;
size_Qs = [num_smp, num_cluster, num_view];
[Qs,F,H,G,C,Es,J,Theta,beta_s,Eij] = Initi(in_X,num_view,k);
tempF = F;
disp('Updateing......')
for iter = 1 : max_iter
    S = update_S(F,Qs,C,Es,Theta,beta_s,delta,Eij);
    C = update_C(S,Es,Theta,H,G,num_smp,num_view,alpha,delta);
    Es = update_E(S,C,Theta,delta);
    tempH = alpha * C * G + G;
    [Ut, ~, Vt] = svds(tempH , k);
    H = Ut * Vt';
    tempG = alpha * C' * H + H;
    [Ug, ~, Vg] = svds(tempG , k);
    G = Ug * Vg';
    for idx = 1 : num_view
        St = S{idx};
        Qt = Qs{idx};
        temp_3 = St * Qt;
        beta = beta_s(idx);
        tempF = tempF + beta * temp_3;
    end
    [Uf, ~, Vf] = svds(tempF , k);
    F = Uf * Vf';
    F = max(F,0);
    Qs = update_Q(S,F,J,beta_s);
    Q_tensor = cat(num_view, Qs{:,:});
    q = Q_tensor(:);
    [j, ~] = wshrinkObj_weight(q, lambda, size_Qs, 1, 3);
    J_tensor = reshape(j, size_Qs); 
    for i = 1 : num_view
        J{i} = J_tensor(:, :, i);
    end
    for idx = 1 : num_view
        temp_SCE = S{idx} - C - Es{idx};
        temp_Theta = delta * temp_SCE;
        Theta{idx} = Theta{idx} + temp_Theta;
    end
    beta_s = update_beta(S,F,Qs);
    U = [H,F];
end
U = real(U);
Y = kmeans(U,num_cluster,'maxiter',KMAXiter,'replicates',KREPlic,'EmptyAction','singleton');
res = zeros(1, 8);
res(1, :) = Clustering8Measure(true_labs, Y);
acc = res(1, 7);
nmi = res(1, 4);
ari = res(1, 5);
fprintf('%.4f\n',acc);
fprintf('%.4f\n',ari);
fprintf('%.4f\n',nmi);
end