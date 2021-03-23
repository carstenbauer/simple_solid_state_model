% Calculates the following 1D thermodynamic properties depending on the
% used potential:
% Energy E, Potential Energy V, Kinetic Energy T_k, Pressure P

function [E, V, T_k, P] = thermos(x,v,potential,params)

T_k = 1/2*params.m*sum(v.^2,2);

if strcmpi(potential, 'inverse')
    pos_diff = diff(x(:,1:params.N)');
    V = sum(params.k ./ pos_diff).';
    P = (-1)*params.k*1./(x(:,2)-x(:,1))+params.a; 
elseif strcmpi(potential, 'harm')
    V = sum(1/2*params.k*(x(:,2:params.N) - x(:, 1:(params.N-1))).^2,2);
    P = (-1)*params.k*(x(:,2)-x(:,1))+params.a;
elseif strcmp(potential,'lenardJones')
    T = sum(0.5.*realpow(abs(v),2),2);
    V = sum(-params.A*sum(1./realpow(diff(x,1,2),6),2)+params.B*sum(1./realpow(diff(x,1,2),12),2),2);
    E = T+V;
    P = zeros(size(E));
end

E = T_k+V;

%V=V-V(1)