% Consumes the simulations inital settings and delivers the systems time
% evolution.

function [t u_total] = timeevolution2D(potential, solver, init, params)

% Anonymous functions:
g_harm = @(t, u) harm2D(t, u, params);
g_inverse = @(t, u) inverse2D(t, u, params);
g_comb = @(t, u) comb2D(t, u, params);

%Vectorize along the columns
u0 = state_enc(init.u);

if strcmp(solver,'eulerRich')
    solve = @eulerRich;
else
    solve = @ode45;
end


if (strcmp(potential, 'harm'))
    if strcmp(solver, 'eulerRichAdapt')
        [t, u_total] = eulerRichAdapt(u0, g_harm, params.t_i, params.dt, params.t_final,params.RelTol);
    else
        [t, u_total] = solve(g_harm, params.t_i:params.dt:params.t_final, u0);
    end
elseif (strcmp(potential, 'inverse'))
    if strcmp(solver, 'eulerRichAdapt')
        [t, u_total] = eulerRichAdapt(u0,g_inverse, params.t_i,params.dt,params.t_final, params.RelTol);
    else
        [t, u_total] = solve(g_inverse, params.t_i:params.dt:params.t_final, u0);
    end
elseif (strcmp(potential, 'comb'))
    if strcmp(solver, 'eulerRichAdapt')
        [t, u_total] = eulerRichAdapt(u0,g_comb, params.t_i,params.dt,params.t_final, params.RelTol);
    else
        [t, u_total] = solve(g_comb, params.t_i:params.dt:params.t_final, u0);
    end    
end
