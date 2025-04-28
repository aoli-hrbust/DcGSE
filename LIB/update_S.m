function Ss = update_S(F,Qs,C,Es,Theta,beta_s,delta,Eij)
num_smp = size(F,1);
num_view = length(Qs);
Ss = cell(1,num_view);
temp_Sij = zeros(num_smp,num_smp);
for v = 1 : num_view
    temp_Qsv = Qs{v};
    temp_Es = Es{v};
    temp_Theta = Theta{v};
    beta_v = beta_s(v);
    temp_FQv = F * temp_Qsv';
    temp_S3 = 1 + beta_v + delta;
    temp_e = Eij{v};
    for i = 1 : num_smp
        for j = 1 : num_smp
            temp_Cij = C(i,j);
            temp_Eij = temp_Es(i,j);
            temp_Tij = temp_Theta(i,j);
            temp_fqvij= temp_FQv(i,j);
            temp_eij = temp_e(i,j);
            temp_S1 = delta * (temp_Cij + temp_Eij);
            temp_S2 = (beta_v *  temp_fqvij + temp_S1 - temp_Tij - 0.5 * temp_eij);
            temp_sij = temp_S2 / temp_S3;
            temp_Sij(i,j) = temp_sij;
        end
    end
    Sij = (temp_Sij + temp_Sij')/2;
    Sij = Sij - diag(diag(Sij));
    Ss{v} = Sij;
end