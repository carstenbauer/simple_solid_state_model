clc;

%% -- Constants
params.N = 20;          % number of particles
params.a = 1;           % initial lattice spacing
params.k = 100;         % spring constant
params.m = 1;           % mass
V = 1;                  % maximal initial velocities amplitude
                        % (velocities in [-V,V])
                        
params.t_final = 100;   % simulation end time
params.dt = 0.01;       % stepsize in time

potential = 'harm';  % used potential: 'harm', 'inverse', 'lenardJones'
params.A = 2;           % lennard-jones attractive
params.B = 64;          % lennard-jones repulsive

solver = 'ode45';       % solver to use: 'eulerRich', 'eulerRichAdapt', 'ode45'
params.RelTol = 1E-3;   % in case of eulerRichAdapt

params.bins = 100;      % nbins: degree of discretization
params.avrg_period = 2; % average period for temperature calculation

% params.piston = 1;   % piston trigger
% params.DeltaL = 10;  % penetration depth
% params.tau = 40;     % penetration time

params.oscWall = 1;   % oscillating wall trigger


% Internal
TitleFontSize = 12;



%% -- Init

 % inital particle positions
init.x = ((1:params.N)-1)*params.a;

% inital random velocity distribution with borders fixed (v=0)
init.v = [0 random('unif',-V,V,1,params.N-2) 0];
                                                   
disp('Init complete');
%%

from = 0.1;
step = 0.1;
till = 9;

w=from:step:till;
y=zeros(size(w));

for k = 1:length(w)
    params.omega = w(k);

%% -- Solve

% ode45
tic
[t,x,v] = timeevolution(potential,solver,init,params);
dur_ode = toc;
disp(['Solving omega=', num2str(params.omega)]);

%% -- Thermodynamic Properties

[E,V,T_k,P] = thermos(x,v,potential,params);

y(k)=max(E);


%% -- Visualize


% % Static Time Evolution
% figure;
% plot(x,t);
% title('Time-Evolution of Particles Absolute Positions','FontSize',TitleFontSize);
% ylabel('Time t');
% xlabel('Position x');
% ylim([0,params.t_final]);

% % Energies vs time
% figure;
% hold on;
% plot(t,T_k,'c','LineWidth',2);
% plot(t,V,'b','LineWidth',2);
% plot(t,E,'k','LineWidth',2);
% title('Time Evolution of Energies','FontSize',TitleFontSize);
% xlabel('Time');
% ylabel('Energy');
% legend('T_k', 'V', 'E');
% xlim([0,params.t_final]);
% hold off;
end