% --------------------------------------------------------------------
% Basic script for simulating the 2D solid state with focus on motion and
% energy!
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
params.N_x = 20;            % particle count in x direction
params.N_y = 20;            % particle count in y direction
params.t_i = 0;             % simulation starting time
params.t_final = 50;        % simulation final time
params.dt = 0.01;           % simulation stepsize

params.movie_speed = 0.01;  % this scalar specifies the movie speed
                            % (bigger means faster)
                            
potential = 'harm';         % used potential: 'harm', 'inverse', 'comb'
params.sec_neighbors = 0;   % enable more than neighbor interaction
                            % (only for spring potential)

solver = 'ode45';           % solver to use
%params.RelTol = 1E-3;      % in case of eulerRichAdapt


% Internal
TitleFontSize = 12;

%% -- Init

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
disp(['Solving: ', num2str(dur_ode), 's']);


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

% Motion Movie
h=figure;
movegui(h);
jstep = round(size(t,1)/(1/(params.movie_speed)));
z=0;
for j=1:jstep:length(t)
    
    temp_mat = state_dec(u_total(j,:),params.N_y);
    temp_mat = temp_mat(:,1:params.N_x);
    scatter(real(temp_mat(:)), imag(temp_mat(:)),25,[0.8,0,0],'filled');
    %plot(temp_mat, '.','MarkerSize',15);
    
    ylim([-1 (params.N_y+2)]);
    xlim([-1 params.N_x+2]);
    %title('Particles slipping away 2D');
    %axis equal;
    
%     filename = sprintf('image_%0.4d.png',z);
%     saveas(gcf,filename,'png');
    z=z+1;
    
    M(j) = getframe;
end