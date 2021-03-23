% This script was used to create the comparison plot between the harmonic
% spring potential implemented using a for loop and the vectorized approach


clear;
close all;
clc;

%% -- Constants
params.a = 1;           % initial lattice spacing in x-direction
params.k = 1;           % spring constant
params.m = 1;           % mass
V = 0.1;                % order of magnitude of the initial velocities
params.t_final = 50;    % simulation end time
params.dt = 0.01;        % stepsize in time

potential = 'harm';     % used potential
params.a = 2; %lennard-jones attractive
params.b = 64; %lennard-jones repulsive

solver = 'ode45';
params.RelTol = 1E-3; %eulerRichAdapt

% Internal
TitleFontSize = 12;

from=10;
step=100;
till=1910;
odet = zeros(length(from:step:till),1);
odeFORt = zeros(length(from:step:till),1);
 z=1;
for k = from:step:till
    %% -- Init
    params.N=k;           % number of particles in x-direction (horizontal)
    
     % inital particle positions
    init.x = ((1:params.N)-1)*params.a;

    % inital random velocity distribution with borders fixed (v=0)
    init.v = [0, 10*V*rand(1,params.N-2).*sign(randn(1,params.N-2)) , 0];

    %% -- Solve
    % ode45
    tic
    [t_ode,x_ode,v_ode] = timeevolution(potential,'ode45',init,params);
    odet(z) = toc;


    df = @(t, y) harmFOR(y,t,params);
    tic
    [t,x,v] = ode45(df,0:params.dt:params.t_final,[init.x, init.v]);
    odeFORt(z) = toc;

    disp(z);
    z=z+1;
end

%%
%Comparison
figure;
plot(from:step:till,odet,'r');
hold on;
plot(from:step:till,odeFORt,'k');
legend('Vectorized','For loop');
xlabel('Solving Time');
ylabel('Particle Count');