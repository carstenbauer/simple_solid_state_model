% Consumes the simulations inital settings and delivers the systems time
% evolution. Thereby rows correspond to timesteps and columns to particles

function [t, x, v] = timeevolution(potential,solver,init,params)
y_init = [init.x, init.v];

if strcmp(potential,'inverse')
    df = @(t, y) inverse(y,t,params);
elseif strcmp(potential,'lenardJones')
    df = @(t, y) lenardJones(y,t,params);
else
    df = @(t, y) harm(y,t,params);
end


if strcmp(solver, 'eulerRich')
   [t, y] = eulerRich(df,0:params.dt:params.t_final,y_init); 
elseif strcmp(solver, 'eulerRichAdapt')
   [t, y] = eulerRichAdapt(y_init,df,0,params.dt,params.t_final,params.RelTol);
else
   [t, y] = ode45(df,0:params.dt:params.t_final,y_init); 
end

x = y(:,1:params.N);
v = y(:,params.N+1:end);