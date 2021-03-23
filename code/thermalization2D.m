% --------------------------------------------------------------------
% Thermalization script for simulating the 2D solid state with focus on its
% temperature trend (thermalization)
% --------------------------------------------------------------------

clear;
close all;
clc;

%% -- Constants

params.k_spring = 1;        % spring constant harm potential
params.k_coulomb = -1.5;    % spring constant coulomb potential
params.m = 1;               % mass
params.v_amp = 1;           % maximal initial velocity amplitude
                            % velocities in [-v_amp,v_amp]
params.N_x = 30;            % particle count in x direction
params.N_y = 30;            % particle count in y direction
params.t_i = 0;             % simulation starting time
params.t_final = 100;       % simulation final time
params.dt = 0.1;            % simulation stepsize

params.movie_speed = 0.01;  % this scalar specifies the movie speed
                            % (bigger means faster)
                            
potential = 'harm';         % used potential: 'harm', 'inverse', 'comb'
params.sec_neighbors = 0;   % enable more than neighbor interaction
                            % (only for spring potential)

solver = 'ode45';           % solver to use
%params.RelTol = 1E-3;      % in case of eulerRichAdapt


params.bins = 100;      % nbins: degree of discretization
params.avrg_period = 2; % average period for temperature calculation

% Internal
TitleFontSize = 12;

%% -- Init
% Note: u0 can be either a row or a column vector
% Use complex numbers to represent two directions
init.u = zeros(params.N_y,2*params.N_x);
% Assign the initial positions to the matrix u0 in the format z= x + i * y
pos = repmat((1:params.N_x), [params.N_y 1]) + (repmat((i*(1:params.N_y)),[params.N_x 1])).';
vel = params.v_amp*(rand((params.N_y-2),(params.N_x-2))-0.5 + i*(rand((params.N_y-2),(params.N_x-2))-0.5));
%vel = params.v_amp * (rand((params.N_y-2), (params.N_x-2))-0.5); % only in x-direction

%vel = zeros((params.N_y-2), (params.N_x-2));
%vel(1, round(params.N_x/2)) = i*params.v_amp;
%vel(1, round(params.N_x/2)-1) = i*params.v_amp;

init.u(:, 1:params.N_x) = pos;
init.u(2:(params.N_y-1), (params.N_x+2):(2*params.N_x-1)) = vel;

disp('Init complete');

%% -- Solve
disp('Solving...');

tic;
[t, u_total] = timeevolution2D(potential, solver, init, params);
dur_ode = toc;
disp(['Ode45 solving: ', num2str(dur_ode), 's']);


%% -- Thermodynamic Properties

[E, V, T_k, P] = thermos2D(u_total, potential, params);

disp('Energies calculated');


%% -- Visualize

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
xlim([params.t_i,params.t_final]);
hold off;


e_i = zeros(size(u_total,1),params.N_x*params.N_y);
for k = 1:size(u_total,1)
    mat_temp =  state_dec(u_total(k,:), params.N_y);
    v = abs(mat_temp(:, (params.N_x + 1):end));
    e_i_t = 1/2*params.m*v.^2;
    e_i_t = e_i_t(:).'; % including borders with v=0
    e_i(k,:) = e_i_t;
end
%e_i now stores the kinetic energy per particle (columns) at each time
%(row)


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
    [T, energy_bins, n] = temperature(e_i(times(k),:),params);
    
    subplot(2,2,k), bar(energy_bins,n,'b','EdgeColor','none');
    subplot(2,2,k), title(['t=' num2str(t(times(k))) 's']);
    subplot(2,2,k), alpha(0.3);
    
    theo = exprnd(T,1,params.N_x*params.N_y);
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


% Motion Movie
h=figure;
movegui(h);
jstep = round(size(t,1)/(1/(params.movie_speed)));
for j=1:jstep:length(t)
    
    %Alternative plots:
    temp_mat = state_dec(u_total(j,:),params.N_y);
    temp_mat = temp_mat(:,1:params.N_x);
    %scatter(real(temp_mat(:)), imag(temp_mat(:)),25,[0.8,0,0],'filled');
    plot(temp_mat, '.','MarkerSize',15);
    
    ylim([-1 (params.N_y+2)]);
    xlim([-1 params.N_x+2]);
    %title('Movement of particles on a rectangle');
    %axis equal;
    M(j) = getframe;
end



