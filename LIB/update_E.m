function Es = update_E(S,C,Theta,delta)
num_view = length(S);
Es = cell(1,num_view);
tau = 1/delta;
num_view = length(S);
for idx = 1 : num_view
    temp = S{idx} - C + Theta{idx}/delta;
    temp_Es = soft(temp, tau);
    Es{idx} = temp_Es;
end
end