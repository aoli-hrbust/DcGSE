<<<<<<< HEAD
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
=======
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
>>>>>>> e1248d18d20a74f08f96760f80be9bf1df35eab1
end