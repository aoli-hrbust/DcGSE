function beta_s = update_beta(S,F,Qs)
num_view = length(S);
beta_s = zeros(num_view, 1);
nu = zeros(num_view, 1);
x = 0;
for v = 1:num_view
    Mtmp = S{v} - F * Qs{v}';
    nu(v) = trace(Mtmp * Mtmp'); %vp
    x = x + 1 / nu(v); 
end
for v = 1:num_view
    beta_s(v) = 1 / (x * nu(v)) ^ 2;
end
end
