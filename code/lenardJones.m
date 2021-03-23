% Implementation of the derivative of the lennard-jones potential
% (supporting piston + oscWall)

function dy = lenardJones(y0,t,params)

dy(1) = 0; %derivative position left boundary
dy(params.N) = 0; %derivative position right boundary
dy(params.N+1) = 0; %derivative velocity left boundary
dy(2*params.N) = 0; %derivative velocity right boundary

if (isfield(params,'piston')) && (params.piston) && (t<params.tau)
    dy(1) = 1/params.tau*params.DeltaL;
end

if (isfield(params,'oscWall')) && (params.oscWall)
    dy(1) = cos(params.omega*t);
end

dy(2:(params.N-1)) = y0((params.N+2):(end-1));

for n = (params.N+2):(2*params.N-1)
    m = n-params.N;
    d1 = y0(m-1)-y0(m);
    d2 = y0(m+1)-y0(m);
    dy1 = -6*params.A*(1/abs(d1))^7 + 12*params.B*(1/abs(d1))^13;
    dy2 = +6*params.A*(1/abs(d2))^7 - 12*params.B*(1/abs(d2))^13;
    dy(n) = dy1 + dy2; 
end

dy = dy.'; %return column vector