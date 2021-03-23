% This bunchy script was for measuring the critical velocities for
% different potentials and solvers in the sense that particles slip away


clear;
close all;
clc;

%% -- Constants
params.N = 20;           % number of particles in x-direction (horizontal)
params.a = 1;           % initial lattice spacing in x-direction
params.k = 1;           % spring constant
params.m = 1;           % mass

params.t_final = 10;    % simulation end time
params.dt = 0.01;        % stepsize in time

potential = 'harm';     % used potential

params.RelTol = 1E-3; %eulerRichAdapt

from = 5;
step = 0.1;
till = 10;

%first harm second inverse
odeverr = zeros(5,2);
% erverr = zeros(5,2);
eraverr = zeros(5,2);

for o = 1:5

for l = 2:2
    
    if l == 1
        potential = 'harm';
    elseif l ==2
        potential = 'inverse';
    end

for V = from:step:till

     % inital particle positions
    init.x = ((1:params.N)-1)*params.a;

    % inital random velocity distribution with borders fixed (v=0)
    init.v = [0 random('unif',-V,V,1,params.N-2) 0];

    [t,x_ode,v] = timeevolution(potential,'ode45',init,params);
%     [t,x_er,v] = timeevolution(potential,'eulerRich',init,params);
    [t,x_era,v] = timeevolution(potential,'eulerRichAdapt',init,params);

    
    if (max(max(x_ode))>init.x(params.N) || min(min(x_ode))<init.x(1)) && odeverr(o,l)== 0
        odeverr(o,l) = V;
    end
    
%     if ((max(max(x_er))>init.x(params.N) || min(min(x_er))<init.x(1))) && erverr(o,l)== 0
%         erverr(o,l) = V;
%     end
    
    if ((max(max(x_era))>init.x(params.N) || min(min(x_era))<init.x(1))) && eraverr(o,l)== 0
        eraverr(o,l) = V;
    end

end
end
end