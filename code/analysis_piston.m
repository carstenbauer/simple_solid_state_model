% This script was used to analyse the systems reaction on a piston. First
% a thermalized system is created. Afterwards the piston simulation will be
% done. Focus in plots lies on energy and temperature

clear;
close all;
clc;

%% -- Constants
params.N = 1000;           % number of particles in x-direction (horizontal)
params.a = 1;           % initial lattice spacing in x-direction
params.k = 1;           % spring constant
params.m = 1;           % mass
V = 1;                  % order of magnitude of the initial velocities
params.t_final = 200;    % simulation end time
params.t_final_piston = 400; 
params.dt = 0.1;        % stepsize in time

solver = 'ode45';        % used solver

potential = 'harm';      % used potential

params.bins = 100;      % bincount
params.avrg_period = 2; % average period for Temperature

params.DeltaL = 300;
params.tau = 280;

% Internal
TitleFontSize = 12;


%% -- Init
init.x = ((1:params.N)-1)*params.a;
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


% % % Energies vs time
% % figure;
% % hold on;
% % plot(t,T_k,'c','LineWidth',2);
% % plot(t,V,'b','LineWidth',2);
% % plot(t,E,'k','LineWidth',2);
% % title('Time Evolution of Energies','FontSize',TitleFontSize);
% % xlabel('Time');
% % ylabel('Energy');
% % legend('T_k', 'V', 'E');
% % xlim([0,params.t_final]);
% % hold off;
% 
% 
% % Kinetic Energy per particle matrix (row = timestep)
% e_i = 1/2*params.m*v.^2;
% 
% 
% disp('Calculating Temperature Evolution...');
% tic;
% % Temperature Time Evolution
% [T, energy_bins, n, quality] = temperature(e_i, params);
% count = round(params.avrg_period / params.dt);
% t_plot = t(1:count:end);
% T_plot = zeros(size(t_plot));
% quality_plot = zeros(size(t_plot));
% idx = 1;
% for k = count:count:length(t)
%     T_plot(idx) = mean(T(k-count+1:k));
%     quality_plot(idx) = mean(quality(k-count+1:k));
%     t_plot(idx) = t(k);
%     idx=idx+1;
% end
% t_plot = [t(1); t_plot(1:end-1)];
% T_plot = [T(1); T_plot(1:end-1)];
% quality_plot = [quality(1); quality_plot(1:end-1)];
% 
% figure;
% subplot(1,2,1), plot(t_plot,T_plot,'b');
% xlabel('Time');
% ylabel('Temperature');
% subplot(1,2,2), plot(t_plot,quality_plot,'b');
% xlabel('Time');
% ylabel('Temperature Quality');
temp = toc;
disp(['Done. (', num2str(temp) ,'s)']);

% ---------------------- THERMALIZED -------------------------------
%%
init.v = v(end,:);
init.x = x(end,:);

params.piston = 1;
params.t_final = params.t_final_piston;

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
% times = [1, round(length(t)/3), round(length(t)*2/3), length(t)];
% 
% figure;
% for k = 1:4
%     [T, energy_bins, n, quality] = temperature(e_i(times(k),:),params);
%     
%     subplot(2,2,k), bar(energy_bins,n,'b','EdgeColor','none');
%     subplot(2,2,k), title(['t=' num2str(t(times(k))) 's']);
%     subplot(2,2,k), alpha(0.3);
%     
%     theo = exprnd(T,1,params.N);
%     [n_theo energy_bins] = hist(theo,energy_bins);
% 
% %     n_theo = params.N*Z^(-1)*exp(-energy_bins/T);
% 
%     subplot(2,2,k), hold on;
%     subplot(2,2,k), bar(energy_bins, n_theo,'r','EdgeColor','none');
%     subplot(2,2,k), alpha(0.3);
% 
%     
%     
%     if k == 4
%         subplot(2,2,k), title(['t=' num2str(t(times(k))) 's, T=', num2str(T)]);
%     end
% end

temp = toc;
disp(['Done. (', num2str(temp) ,'s)']);