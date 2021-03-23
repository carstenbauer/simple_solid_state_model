% --------------------------------------------------------------------
% Boundaries script for simulating the 1D solid states behavior under
% external stimulation
% --------------------------------------------------------------------

clear;
close all;
clc;

%% -- Constants
params.N = 20;          % number of particles
params.a = 1;           % initial lattice spacing
params.k = 10;         % spring constant
params.m = 1;           % mass
V = 1;                  % maximal initial velocities amplitude
                        % (velocities in [-V,V])
                        
params.t_final = 100;   % simulation end time
params.dt = 0.01;       % stepsize in time

potential = 'inverse';  % used potential: 'harm', 'inverse', 'lenardJones'
params.A = 2;           % lennard-jones attractive
params.B = 64;          % lennard-jones repulsive

solver = 'ode45';       % solver to use: 'eulerRich', 'eulerRichAdapt', 'ode45'
params.RelTol = 1E-3;   % in case of eulerRichAdapt

params.bins = 100;      % nbins: degree of discretization
params.avrg_period = 2; % average period for temperature calculation

params.piston = 1;   % piston trigger
params.DeltaL = 10;  % penetration depth
params.tau = 70;     % penetration time

% params.oscWall = 1;   % oscillating wall trigger
% params.omega = 0.5;     % oscillation frequency

% Internal
TitleFontSize = 12;

%% -- Init

 % inital particle positions
init.x = ((1:params.N)-1)*params.a;

% inital random velocity distribution with borders fixed (v=0)
init.v = [0 random('unif',-V,V,1,params.N-2) 0];
                                                   
disp('Init complete');
%% -- Solve
disp('Solving...');

% ode45
tic
[t,x,v] = timeevolution(potential,solver,init,params);
dur_ode = toc;
disp(['Solving: ', num2str(dur_ode), 's']);

%% -- Thermodynamic Properties

[E,V,T_k,P] = thermos(x,v,potential,params);
disp('Energies calculated');


%% -- Visualize


% Static Time Evolution
figure;
plot(x,t);
title('Time-Evolution of Particles Absolute Positions','FontSize',TitleFontSize);
ylabel('Time t');
xlabel('Position x');
ylim([0,params.t_final]);

% Energies vs time
figure;
hold on;
plot(t,T_k,'c','LineWidth',2);
plot(t,V,'b','LineWidth',2);
plot(t,E,'k','LineWidth',2);
title('Time Evolution of Energies','FontSize',TitleFontSize);
xlabel('Time');
ylabel('Energy');
legend('T_k', 'V', 'E');
xlim([0,params.t_final]);
hold off;

% Sub-Plot Animation
h = figure;
movegui(h); %% DUAL Monitor BUG
jstep = round(size(t,1)/200);
z=0;
for j=1:jstep:size(t)
    subplot(2,1,1), plot(x,t);
    subplot(2,1,1), hold on;
    subplot(2,1,1), plot(x(j,:),t(j),'black.');
    title('Animation of the Systems Time Evolution','FontSize',TitleFontSize);
    ylabel('Time t');
    xlabel('Position x');
    xlim([(min(min(init.x))-2*params.a),(max(max(init.x))+2*params.a)]);
    ylim([min(min(t)), max(max(t))]);
    
    hold off;
    subplot(2,1,2), plot(x(j,:),zeros(params.N),'blacko','LineWidth',2,'MarkerFaceColor','k');
    xlabel('Position x');
    xlim([(min(min(x))-2*params.a),(max(max(x))+2*params.a)]);
    set(gca, 'ytick', []);
    
%     filename = sprintf('image_%0.4d.png',z);
%     saveas(gcf,filename,'png');
%     z=z+1;
    
    M(j) = getframe;
end
