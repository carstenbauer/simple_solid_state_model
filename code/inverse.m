% Implementation of the derivative of the inverse coulomb potential
% (supporting piston + oscWall)

function dF = inverse(y,t,params)

% Boundaries
dF(1) = 0; %derivative position left boundary
dF(params.N) = 0; %derivative position right boundary
dF(params.N+1) = 0; %derivative velocity left boundary e.g. external force
dF(2*params.N) = 0; %derivative velocity right boundary

if (isfield(params,'piston')) && (params.piston) && (t<params.tau)
    dF(1) = 1/params.tau*params.DeltaL;
end

if (isfield(params,'oscWall')) && (params.oscWall)
    dF(1) = cos(params.omega*t);
end


dF(2:(params.N-1)) = y((params.N+2):(end-1));

for n = (0.5*length(y)+2):(length(y)-1)
    m = n-0.5*length(y);
    dF(n) = params.k*(-(1/(y(m)-y(m+1))).^2 + (1/(y(m)-y(m-1))).^2); 
end

dF = dF.'; %return column vector