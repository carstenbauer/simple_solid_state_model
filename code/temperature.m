% Given the kinetic energies per particle (can be a matrix where one row
% corresponds to one timestep) the function calculates the temperature by
% fitting the cumulative distribution function. Moreover it returns the
% discretized energies (energy_bins), the particle count per bin (n) and a
% temperature quality measure (quality - usually in between 1 and 10 - less
% means good)

function [T, energy_bins, n, quality] = temperature(e_i,params)

T=zeros(size(e_i,1),1);
quality=zeros(size(e_i,1),1);
energy_bins = zeros(size(e_i,1),params.bins);

for k = 1:size(e_i,1)
    [n(k,:) energy_bins(k,:)] = hist(e_i(k,:), params.bins);
    n_cumul = zeros(1, length(n(k,:)));
    for u = 1:length(n(k,:))
        n_cumul(u) = sum(n(k,1:u));
    end
    
    if isfield(params,'N')
        n_cumul_norm = n_cumul / params.N;
    else
        n_cumul_norm = n_cumul / (params.N_x*params.N_y);
    end

    bestcoeffs=fminsearch(@boltzmann,1,[],energy_bins(k,:),n_cumul_norm);
    lambda = bestcoeffs(1);
    quality(k) = boltzmann(lambda,energy_bins(k,:),n(k,:)) * realpow(10,-5);
    T(k) = 1/lambda;
end