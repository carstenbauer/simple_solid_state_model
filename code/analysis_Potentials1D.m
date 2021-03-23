% Comparison of the motions due to different potentials

clear;
close all;
clc;

%% -- Constants
params.N = 20;           % number of particles in x-direction (horizontal)
params.a = 1;           % initial lattice spacing in x-direction
params.k = 1;           % spring constant
params.m = 1;           % mass
V = 0.1;                % order of magnitude of the initial velocities
params.t_final = 20;    % simulation end time
params.dt = 0.01;        % stepsize in time

params.A = 2; %lennard-jones attractive
params.B = 64; %lennard-jones repulsive

% Internal
TitleFontSize = 12;

%% -- Init

 % inital particle positions
init.x = ((1:params.N)-1)*params.a;

% inital random velocity distribution with borders fixed (v=0)
init.v = [0, 10*V*rand(1,params.N-2).*sign(randn(1,params.N-2)) , 0];
                                                   
disp('Init complete');
%% -- Solve

[t,x_harm,v_harm] = timeevolution('harm','ode45',init,params);
disp('harm done');
[t,x_inv,v_inv] = timeevolution('inverse','ode45',init,params);
disp('coul done');
[t,x_lj,v_lj] = timeevolution('lenardJones','ode45',init,params);
disp('lenjon done');

%% -- Visualize

% Static Time Evolution
h=figure;
movegui(h);

xmin=-5;
xmax=23;

hold on;
subplot(1,3,1), plot(x_harm,t);
xlim([xmin xmax]);
xlabel('Position');
ylabel('Time');
title('Spring');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot(1,3,2), plot(x_inv,t);
xlim([xmin xmax]);
xlabel('Position');
ylabel('Time');
title('Coulomb');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
subplot(1,3,3), plot(x_lj,t);
xlim([xmin xmax]);
xlabel('Position');
ylabel('Time');
title('Lennard-Jones');
set(gca,'xtick',[]);
set(gca,'ytick',[]);
