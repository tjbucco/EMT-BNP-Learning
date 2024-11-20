matfiles = dir("*.mat");
file_inc_prime = 0;
for file = matfiles'
    file_inc_prime = file_inc_prime + 1
    load(file.name);
    %parameters for OSPA metric
    p = 1; cutoff = 120; cutoff_GWD = 1000; cutoff_mse = cutoff*100.*[Elsq_data.^2; Essq_data.^2; ElambdaM_data]; alpha = 2; %OSPA metric parameters

%     meanT(file_inc_prime) = mean(data.Ts - data.Te);
%     meanD(:,:,file_inc_prime) = mean(data.Dstar, 3);
%     meanSig(:,:,file_inc_prime) = mean(data.Sigmavstar, 3);
%     meanLamb(file_inc_prime) = mean(data.lambdaMstar);
%     CG_prime(file_inc_prime) = OSPA_NC_nolabels(end)/OSPA_DP_nolabels(end);

%     [OSPA_DP_newcutoff(file_inc_prime), OSPA_DP_newcutoff_pert(:, file_inc_prime), ~, ~] = calculate_OSPA_notracklabels(data.X, DPresults.X, p, cutoff, alpha);
[~, OSPA_DP_GWD_nolabels_pert_new(:, file_inc_prime), ~, ~, MSE_DP_nolabels_new(:,file_inc_prime), ~, ~]  = calculate_OSPA_notracklabels(data.X, DPresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), DPresults.D, data.q, DPresults.q, data.lambdaMstar(data.C), DPresults.lambdaM, data.C, DPresults.object_class_est.Cest);
% [OSPA_DP(sim), OSPA_DP_pert(:, sim),  ~, ~] = calculate_OSPA(data.X, DPresults.X, data.Ts, data.Te, DPresults.Ts, DPresults.Te, p, cutoff, alpha);
%
%     [OSPA_DP(file_inc_prime), ~,  ~, ~] = calculate_OSPA(data.X, DPresults.X, data.Ts, data.Te, DPresults.Ts, DPresults.Te, p, cutoff, alpha);
%     [OSPA_DP_GWD(file_inc_prime), ~, ~, ~, ~, ~, ~]  = calculate_OSPA(data.X, DPresults.X, data.Ts, data.Te, DPresults.Ts, DPresults.Te, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), DPresults.D, data.q, DPresults.q, data.lambdaMstar(data.C), DPresults.lambdaM, data.C, DPresults.object_class_est.Cest);
    %
    %         %% No clustering metrics
%     [OSPA_NC_newcutoff(file_inc_prime), OSPA_NC_newcutoff_pert(:, file_inc_prime), ~, ~] = calculate_OSPA_notracklabels(data.X, NCresults.X, p, cutoff, alpha);
[~, OSPA_NC_GWD_nolabels_pert_new(:, file_inc_prime), ~, ~, MSE_NC_nolabels_new(:,file_inc_prime), ~, ~]  = calculate_OSPA_notracklabels(data.X, NCresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), NCresults.D, data.q, NCresults.q, data.lambdaMstar(data.C), NCresults.lambdaM, data.C, NCresults.object_class_est.Cest);
% [OSPA_NC(sim), OSPA_NC_pert(:, sim),  ~, ~] = calculate_OSPA(data.X, NCresults.X, data.Ts, data.Te, NCresults.Ts, NCresults.Te, p, cutoff, alpha);
%     %
%     [OSPA_NC(file_inc_prime), ~,  ~, ~] = calculate_OSPA(data.X, NCresults.X, data.Ts, data.Te, NCresults.Ts, NCresults.Te, p, cutoff, alpha);
%     [OSPA_NC_GWD(file_inc_prime), ~, ~, ~, ~, ~, ~]  = calculate_OSPA(data.X, NCresults.X, data.Ts, data.Te, NCresults.Ts, NCresults.Te, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), NCresults.D, data.q, NCresults.q, data.lambdaMstar(data.C), NCresults.lambdaM, data.C, NCresults.object_class_est.Cest);

[~, OSPA_PMBM_GWD_nolabels_pert_new(:, file_inc_prime), ~, ~, MSE_PMBM_nolabels_new(:,file_inc_prime), ~, ~]  = calculate_OSPA_notracklabels(data.X, PMBMresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), PMBMresults.D, data.q, repmat(hyperparams.Essq_data,1, PMBMresults.maxNactive), data.lambdaMstar(data.C), repmat(hyperparams.ElambdaM, 1, PMBMresults.maxNactive),data.C, 1:PMBMresults.maxNactive, 1);

[~, OSPA_FUSION_GWD_nolabels_pert_new(:, file_inc_prime), ~, ~, MSE_FUSION_nolabels_new(:,file_inc_prime), ~,~]  = calculate_OSPA_notracklabels(data.X, FUSIONresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), FUSIONresults.D, data.q, FUSIONresults.q, data.lambdaMstar(data.C), repmat(hyperparams.ElambdaM, 1, FUSIONresults.maxNactive), data.C, FUSIONresults.object_class_est.Cest, 1);

[~, OSPA_Perf_GWD_nolabels_pert_new(:, file_inc_prime), ~, ~, MSE_Perf_nolabels_new(:,file_inc_prime), ~, ~]  = calculate_OSPA_notracklabels(data.X, Perfresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), Perfresults.D, data.q, Perfresults.q, data.lambdaMstar(data.C), Perfresults.lambdaM, data.C, Perfresults.object_class_est.Cest);
% [OSPA_Perf(sim), OSPA_Perf_pert(:, sim),  ~, ~] = calculate_OSPA(data.X, Perfresults.X, data.Ts, data.Te, Perfresults.Ts, Perfresults.Te, p, cutoff, alpha);

end
