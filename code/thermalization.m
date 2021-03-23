% --------------------------------------------------------------------
% Thermalization script for simulating the 1D solid states thermal behavior
% --------------------------------------------------------------------

clear;
clc;

%% -- Constants
params.N = 1000;        % number of particles
params.a = 1;           % initial lattice spacing
params.k = 1;           % spring constant
params.m = 1;           % mass
V = 1;                  % maximal initial velocities amplitude
                        % (velocities in [-V,V])
                        
params.t_final = 500;   % simulation end time
params.dt = 0.1;        % stepsize in time

potential = 'harm';     % used potential: 'harm', 'inverse', 'lenardJones'
params.A = 2;           % lennard-jones attractive
params.B = 64;          % lennard-jones repulsive

solver = 'ode45';       % solver to use: 'eulerRich', 'eulerRichAdapt', 'ode45'
params.RelTol = 1E-3;   % in case of eulerRichAdapt

params.bins = 100;      % nbins: degree of discretization
params.avrg_period = 2; % average period for temperature calculation

% params.piston = 1;    % piston trigger
% params.DeltaL = 300;  % penetration depth
% params.tau = 350;     % penetration time

% params.oscWall = 1;   % oscillating wall trigger
% params.omega = 1;     % oscillation frequency

% Internal
TitleFontSize = 12;


%% -- Init

% Inital particle positions
init.x = ((1:params.N)-1)*params.a;

% Inital random velocity distribution with borders fixed (v=0)
init.v = [0 random('unif',-V,V,1,params.N-2) 0];

% Left half particles in motion, right half standing still
% init.v = zeros(1,params.N);
% init.v(2:round(params.N/2)) =
% 10*V*rand(1,round(params.N/2)-1).*sign(randn(1,round(params.N/2)-1));


disp('Init complete');
%% -- Solve
disp('Solving...');

% ode45
tic
[t,x,v] = timeevolution(potential,solver,init,params);
dur_ode = toc;
disp(['Solving: ', num2str(dur_ode), 's']);


% % Static Time Evolution
% figure;
% plot(x,t);
% title('Time-Evolution of Particles Absolute Positions','FontSize',TitleFontSize);
% ylabel('Time t');
% xlabel('Position x');
% ylim([0,params.t_final]);
%% -- Thermodynamic Properties

[E,V,T_k,P] = thermos(x,v,potential,params);
disp('Energies calculated');


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


% Kinetic Energy per particle matrix (row = timestep)
e_i = 1/2*params.m*v.^2;


disp('Calculating Temperature Evolution...');
tic;
% Temperature Time Evolution
[T, energy_bins, n, quality] = temperature(e_i, params);
count = round(params.avrg_period / params.dt);
t_plot = t(1:count:end);
T_plot = zeros(size(t_plot));
quality_plot = zeros(size(t_plot));
idx = 1;
for k = count:count:length(t)
    T_plot(idx) = mean(T(k-count+1:k));
    quality_plot(idx) = mean(quality(k-count+1:k));
    t_plot(idx) = t(k);
    idx=idx+1;
end
t_plot = [t(1); t_plot(1:end-1)];
T_plot = [T(1); T_plot(1:end-1)];
quality_plot = [quality(1); quality_plot(1:end-1)];

figure;
subplot(1,2,1), plot(t_plot,T_plot,'b');
xlabel('Time');
ylabel('Temperature');
subplot(1,2,2), plot(t_plot,quality_plot,'b');
xlabel('Time');
ylabel('Temperature Quality');



% Time Evolution (discrete) 4 Fits
times = [1, round(length(t)/3), round(length(t)*2/3), length(t)];

figure;
for k = 1:4
    [T, energy_bins, n, quality] = temperature(e_i(times(k),:),params);
    
    subplot(2,2,k), bar(energy_bins,n,'b','EdgeColor','none');
    subplot(2,2,k), title(['t=' num2str(t(times(k))) 's']);
    subplot(2,2,k), alpha(0.3);
    
    theo = exprnd(T,1,params.N);
    [n_theo energy_bins] = hist(theo,energy_bins);

%     n_theo = params.N*Z^(-1)*exp(-energy_bins/T);

    subplot(2,2,k), hold on;
    subplot(2,2,k), bar(energy_bins, n_theo,'r','EdgeColor','none');
    subplot(2,2,k), alpha(0.3);

    
    
    if k == 4
        subplot(2,2,k), title(['t=' num2str(t(times(k))) 's, T=', num2str(T)]);
    end
end

temp = toc;
disp(['Done. (', num2str(temp) ,'s)']);




% % Time Evolution (all)
% h=figure;
% movegui(h);
% kstep = round(size(t,1)/100);

% for k = 1:kstep:size(e_i,1)
%     
%     bar(energy_bins,n);
%     xlim([0, 1.5]);
%     ylim([0, 70]);
%     
%     Temp = [Temp; T];
%     M(k) = getframe;
% end