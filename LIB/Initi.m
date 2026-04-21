function [Qs,F,H,G,C,Es,J,Theta,beta_s,Eij] = Initi(in_X,num_view,k)
num_smp = size(in_X,2);
beta_s = ones(num_view,1) / num_view;
Theta = cell(1,3);
J = cell(1,3);
Es = cell(1,num_view);
Qs = cell(1,num_view);
Eij = cell(1,num_view);
disp('Initialing......')
S = ConMvSs(in_X , k);
end