% 2D Implementation of the derivative of the harmonic spring potential
% (supporting oscWall)

function du = harm2D(t, u, params)

du=zeros(params.N_y*2*params.N_x,1);

du(1:(params.N_x*params.N_y)) = u((params.N_x*params.N_y + 1):end);

du = state_dec(du, params.N_y);
u = state_dec(u, params.N_y);

% Calculate forces
F_horizontal = params.k_spring .* diff(diff(u(2:(params.N_y-1),1:params.N_x).')).';
F_vertical = params.k_spring .* diff(diff(u(:,2:(params.N_x-1))));

% More neighbors:
F_diagonal_tot = zeros(params.N_y, params.N_x);

if (isfield(params,'sec_neighbors') && (params.sec_neighbors))
    u_shift_right_down = zeros(params.N_y, params.N_x);
    u_shift_right_up = zeros(params.N_y, params.N_x);
    u_shift_left_down = zeros(params.N_y, params.N_x);
    u_shift_left_up = zeros(params.N_y, params.N_x);
    u_shift_right_down(2:params.N_y, 2:params.N_x) = u(1:(params.N_y-1),1:(params.N_x-1));
    u_shift_right_up(1:(params.N_y-1), 2:params.N_x) = u(2:params.N_y,1:(params.N_x-1));
    u_shift_left_down(2:params.N_y, 1:(params.N_x-1)) = u(1:(params.N_y-1),2:params.N_x);
    u_shift_left_up(1:(params.N_y-1), 1:(params.N_x-1)) = u(2:params.N_y,2:params.N_x);
    
    F_diagonal_tot = params.k_spring * (u_shift_right_down + u_shift_right_up + u_shift_left_down + u_shift_left_up - 4 * u(1:params.N_y, 1:params.N_x));
end

du(2:(params.N_y-1), (params.N_x+2):(2*params.N_x-1)) = 1/params.m * (F_horizontal + F_vertical + F_diagonal_tot(2:(params.N_y-1),2:(params.N_x-1)));

if (isfield(params,'oscWall')) && (params.oscWall)
    du(:,1) = cos(params.omega*t);
end

du = state_enc(du);