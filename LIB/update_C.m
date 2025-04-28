function [C] = update_C(S,Es,Theta,H,G,num_smp,num_view,alpha,delta)
    temp1 = alpha * eye(num_smp) + delta * eye(num_smp);
    temp = 0;
    temp_2 = cell(1,num_view);
    for idx = 1 : num_view
        temp_2{idx} = S{idx} - Es{idx} - Theta{idx}/delta;
        temp = temp + temp_2{idx};
    end
    C = inv(temp1) * (alpha * H * G' + delta * temp);
end

