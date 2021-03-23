% 2D Implementation of the derivative of the combined potential

function du = comb2D(t, u, params)

du = harm2D(t, u, params);
du_coul = inverse2D(t, u, params);
du((params.N_x*params.N_y+1):end) = du((params.N_x*params.N_y+1):end) + du_coul((params.N_x*params.N_y+1):end);
% The change in position is done both by harm2D and inverse2D.
% Use it only once!