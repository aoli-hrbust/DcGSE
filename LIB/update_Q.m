function Qs = update_Q(S,F,J,beta_s)
num_view = length(S);
Qs = cell(1,num_view);
for idx = 1:num_view
    temp_son = beta_s(idx) * S{idx}' * F + J{idx};
    qtemp = temp_son / (1+beta_s(idx));
    Qs{idx} = qtemp;
end
end