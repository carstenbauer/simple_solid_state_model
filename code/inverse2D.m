% 2D Implementation of the derivative of the inverse coulomb potential
% (supporting oscWall)

function du = inverse2D(t, u, params)

du=zeros(params.N_y*2*params.N_x,1);

du(1:(params.N_x*params.N_y)) = u((params.N_x*params.N_y + 1):end);

du = state_dec(du, params.N_y);
u = state_dec(u, params.N_y);

% Calculate next-neighbor forces
F_tot = zeros(params.N_y, params.N_x);
delta_vertical = diff(u(:,1:params.N_x)); %Note: lower element - upper element --> alternativ sign of k_coulomb
F_vertical = params.k_coulomb * delta_vertical ./ (abs(delta_vertical)).^3;
F_tot(1:(params.N_y-1), :) = F_vertical;
F_tot(2:params.N_y,:) = F_tot(2:params.N_y,:) + (-1) * F_vertical;
delta_horizontal = (diff((u(:,1:params.N_x)).')).'; %Note: right element - left element --> alternativ sign of k_coulomb
F_horizontal = params.k_coulomb * delta_horizontal ./ (abs(delta_horizontal)).^3;
F_tot(:,1:(params.N_x-1)) = F_tot(:,1:(params.N_x-1)) + F_horizontal;
F_tot(:,2:params.N_x) = F_tot(:,2:params.N_x) + (-1) * F_horizontal;
du(2:(params.N_y-1), (params.N_x+2):(2*params.N_x-1)) = 1/params.m * (F_tot(2:(params.N_y-1),2:(params.N_x-1)));

if (isfield(params,'oscWall')) && (params.oscWall)
    du(:,1) = cos(params.omega*t);
end

du = state_enc(du);