function [state_est, object_class_est] = track_Objects_w_knownParams(init_tracking, data, hyper, script, mcmc, sim, j)

total_no_of_objects = init_tracking.maxNactive;
object_class_est = init_tracking.object_class_est;

%% instantiate variables
Xm = NaN(4, init_tracking.maxTactive, total_no_of_objects);
C = zeros(total_no_of_objects, mcmc.L);
Dstar = NaN(2, 2, total_no_of_objects, mcmc.L);
Qstar = NaN(4, 4, total_no_of_objects, mcmc.L);
lambdaMstar = NaN(total_no_of_objects, mcmc.L);

ell = NaN;
X = NaN(4, script.maxT, total_no_of_objects, mcmc.L);

%% initialization
%[Xcs, Ysq, C, Ccs, Dstar, D, Sigmavstar, Sigmav, Sigmavstarcs, Sigmavstarcs, lambdaMstar, lambdaM, lambdaMstarcs] = initialize_estimates(state_est, daoed_est, data, hyper, script, mcmc);
[state_est, Ccs, Dstarcs, Qstarcs, lambdaMstarcs] = initialize_estimates(init_tracking, data, hyper, script, mcmc);
Xcs = state_est.Xcs;
X(:,:,:,1) = Xcs;

[object_class_est, C, Dstar, Qstar, lambdaMstar, Ccs, Dstarcs, Qstarcs, lambdaMstarcs] = cluster_objects(init_tracking, state_est.Xcs, object_class_est, hyper, script.factorials, mcmc, 1, C, Dstar, Qstar, lambdaMstar);

%% MCMC iterations
for ell = 2:mcmc.L

    if ~mod(ell, 10) && script.statustextoutput
        disp(['simulation ', num2str(sim), ', j = ', num2str(j), ', ell = ', num2str(ell)])
    end

    % sample new X
    if ~mod(ell, mcmc.Xrsi)
        Xcs = track_objects_PGAS(init_tracking.Z, Xcs, init_tracking.Ts, init_tracking.Te, Ccs, Dstarcs, hyper.F, Qstarcs, hyper.Sigmax1, mcmc.R, init_tracking.maxNactive);
        X(:,:,:,ell) = Xcs;
    end

    % sample new object parameters and class indicators
    if ~mod(ell, mcmc.DPrsi)
        [object_class_est, C, Dstar, Qstar, lambdaMstar, Ccs, Dstarcs, Qstarcs, lambdaMstarcs] = match_objects_w_knownparams(init_tracking, data, Xcs, object_class_est, hyper, script.factorials, mcmc, ell, C, Dstar, Qstar, lambdaMstar);
    end
end

%% MMSE and MAP
[unq_samps,~,ic_unq] = unique(C(:, mcmc.bis + 1:end).', 'rows');
map_i = mode(ic_unq);
object_class_est.Cest = unq_samps(map_i,:).';
object_class_est.Dest = squeeze(nanmean(Dstar(:, :, :, mcmc.bis + 1:end), 4));
object_class_est.Qest = squeeze(nanmean(Qstar(:, :, :, mcmc.bis + 1:end), 4));
object_class_est.lambdaMest = squeeze(nanmean(lambdaMstar(:, mcmc.bis + 1:end), 2));
state_est.Xest = mean(X(:,:,:,mcmc.bis + 1:end), 4);