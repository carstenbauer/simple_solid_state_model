% Implementation of the derivative of the harmonic spring potential
% (supporting piston + oscWall)

function dF = harm(y,t,params)
% Expects y to be a vector with the particle position
% before the particle velocities

dF = zeros(length(y),1);

% Change of the positions
dx(1:params.N) = y((params.N+1):end);

% Change of the velocities
x(1:params.N) = y(1:params.N);
x_diff = diff(x);
x_diff_shift = [0, x_diff(1:end-1)];
dv = params.k/params.m * (x_diff-x_diff_shift);

% Boundaries
dv(1) = 0;
dv(params.N) = 0;

if (isfield(params,'piston')) && (params.piston) && (t<params.tau)
    dx(1) = 1/params.tau*params.DeltaL;
end

if (isfield(params,'oscWall')) && (params.oscWall)
    dx(1) = cos(params.omega*t);
end


dF(1:params.N) = dx;
dF(params.N+1:end) = dv;