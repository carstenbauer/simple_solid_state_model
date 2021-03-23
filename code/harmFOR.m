function dF = harmFOR(y,t,params)
% Expects y to be a vector with the particle position
% before the particle velocities

dF = zeros(length(y),1);

% Change of the positions
dx(1:params.N) = y((params.N+1):end);

% Change of the velocities
for j = 1:params.N
    v = 0; % new velocity value
    
    if((j~=1) && (j~=params.N))
        v = params.k/params.m * (y(j+1) - y(j) - (y(j) - y(j-1)));
    end;
    
    dv(j) = v;    
end;

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