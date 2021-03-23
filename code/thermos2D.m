% Calculates the following 2D thermodynamic properties depending on the
% used potential:
% Energy E, Potential Energy V, Kinetic Energy T_k, Pressure P

function [E, V, T_k, P] = thermos2D(u_total, potential, params)
T_k = zeros(1, length(u_total(:,1)));
V = zeros(1, length(u_total(:,1)));
E = zeros(1, length(u_total(:,1)));
P = zeros(1, length(u_total(:,1)));

for j = 1:length(u_total(:,1))
    V_harm = 0;
    V_coul = 0;
    P(j) = 0;
    mat_temp =  state_dec(u_total(j,:), params.N_y);
    T_k(1,j) = 0.5*params.m*sum(sum(abs(mat_temp(:, (params.N_x + 1):end)) .^ 2));
    if ((strcmp(potential, 'harm')) || (strcmp(potential, 'comb')))
        V_horizontal = sum(sum(params.k_spring / 2 * abs(diff(mat_temp(2:(params.N_y-1),1:params.N_x).').').^2));
        V_vertical = sum(sum(params.k_spring / 2 * (abs(diff(mat_temp(:,2:(params.N_x-1)))) .^2)));
        V_harm = V_horizontal + V_vertical;
        P(j) = params.k_spring.*sum(real(mat_temp(:,2)-mat_temp(:,1))-1,1)./(params.N_x-1);
        
        if (isfield(params,'sec_neighbors') && (params.sec_neighbors))
            mat_shift_right_down = mat_temp(1:(params.N_y-1), 1:(params.N_x-1));
            V_diag = sum(sum(0.5*params.k_spring*((abs(mat_temp(2:params.N_y,2:params.N_x) - mat_shift_right_down)).^2)));
            mat_shift_left_down = mat_temp(1:(params.N_y-1), 2:params.N_x);
            V_diag = V_diag + sum(sum(0.5*params.k_spring*((abs(mat_temp(2:params.N_y,1:(params.N_x-1)) - mat_shift_left_down)).^2)));
            V_harm = V_harm + V_diag;
        end
    end
    if ((strcmp(potential, 'inverse')) || (strcmp(potential, 'comb')))
        V_coul = (-1)*sum(sum(params.k_coulomb ./ abs(diff(mat_temp(:,1:params.N_x),1,1))));
        V_coul = V_coul + (-1) * sum(sum(params.k_coulomb ./ abs(diff(mat_temp(:,1:params.N_x),1,2))));
    end
    V(1,j) = (V_harm + V_coul);
end

E = T_k + V;
