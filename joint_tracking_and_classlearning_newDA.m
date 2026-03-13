function [results, OSPA, OSPA_pert, OSPA_GWD, OSPA_GWD_pert, MSE] = joint_tracking_and_classlearning_newDA(data, hyperparams, script_params, mcmc_params, seq_tracking, sim)

for j = 1:script_params.J

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DP-based tracking and clustering stage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [state_est_mcmc, object_class_est_mcmc] = track_and_cluster_Objects(seq_tracking, data, hyperparams, script_params, mcmc_params, sim, j);

    q = NaN(2, seq_tracking.maxNactive);
    for i = 1:seq_tracking.maxNactive
        q(1,i) = object_class_est_mcmc.Qest(1,1,object_class_est_mcmc.Cest(i));
        q(2,i) = object_class_est_mcmc.Qest(3,3,object_class_est_mcmc.Cest(i));
    end
    [OSPA(j), OSPA_pert(:, j)] = calculate_OSPA(data.X, state_est_mcmc.Xest, data.Ts, data.Te, seq_tracking.Ts, seq_tracking.Te, script_params.p, script_params.cutoff, script_params.alpha);
    [OSPA_GWD(j), OSPA_GWD_pert(:, j), MSE(:,j), ~, ~]  = calculate_OSPA(data.X, state_est_mcmc.Xest, data.Ts, data.Te, seq_tracking.Ts, seq_tracking.Te, script_params.p, script_params.cutoff_GWD, script_params.alpha, script_params.cutoff_mse, data.Dstar(:,:, data.C), object_class_est_mcmc.Dest(:,:,object_class_est_mcmc.Cest), data.q, q, data.lambdaMstar(data.C), object_class_est_mcmc.lambdaMest(object_class_est_mcmc.Cest), data.C, object_class_est_mcmc.Cest);

    %% insert initial state estimates as measurements
    data_aug = data;
%     for i = 1:seq_tracking.maxNactive
%         for t = seq_tracking.Ts(i):seq_tracking.Te(i)
%             measurements_tent = data_aug.measurements_all{t};
%             data_aug.measurements_all{t} = cat(2, measurements_tent, state_est_mcmc.Xest(1:2,t,i));
%         end
%     end
    [~, no_of_measurements_cell] = cellfun(@size, data_aug.measurements_all);
    max_no_of_meas = max(no_of_measurements_cell);

    %% calculate empircal average of extent for DAOD stage
    D0_emp = mean(object_class_est_mcmc.Dest(:,:,object_class_est_mcmc.Cest),3);
    hyper_emp = hyperparams;
    hyper_emp.D0 = D0_emp;

    if j < script_params.J
        object_class_est = object_class_est_mcmc;
        Z = NaN(2, size(data.Z, 2), script_params.maxT, script_params.binoN*script_params.maxT + 1);
        exists_and_detected_ind = [];
        for t = 1:script_params.maxT

            if t==1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% initialization
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                state_est = struct('Xest', cat(1, squeeze(state_est_mcmc.Xest(1:2,t,:)), squeeze(state_est_mcmc.Xest(3:4,t,:))./1e3), 'Sigmax', repmat(script_params.SigmaXinit, 1, 1, size(state_est_mcmc.Xest, 3)), 'rotangle', zeros(1, size(state_est_mcmc.Xest, 3)), 'Sigmarot', repmat(script_params.Sigmarot0, 1, 1, size(state_est_mcmc.Xest, 3)), 'driving_noise_est', zeros(4, size(state_est_mcmc.Xest, 3)) ,'no_of_objectstminus1', 0);
                %state_est = struct('Xest', zeros(4, daoed_est.no_of_objectst), 'Sigmax', repmat(SigmaXinit, 1, 1, daoed_est.no_of_objectst), 'rotangle', zeros(1,daoed_est.no_of_objectst), 'Sigmarot', repmat(Sigmarot0, 1, 1, daoed_est.no_of_objectst), 'driving_noise_est', zeros(4, daoed_est.no_of_objectst) ,'no_of_objectstminus1', 0);
                object_class_est = struct('Dest', object_class_est_mcmc.Dest(:,:,object_class_est_mcmc.Cest), 'Qest', object_class_est_mcmc.Qest(:,:,object_class_est_mcmc.Cest), 'Cest', object_class_est_mcmc.Cest, 'lambdaMest', object_class_est_mcmc.lambdaMest(object_class_est_mcmc.Cest), 'existsteps', zeros(1, size(state_est_mcmc.Xest, 3)), 'Msum', zeros(1, size(state_est_mcmc.Xest, 3)), 'Zsqsum', zeros(2, size(state_est_mcmc.Xest, 3)), 'Vsqsum', zeros(4, size(state_est_mcmc.Xest, 3)));

                daoed_est.exist_ind_all = zeros(1, size(state_est_mcmc.Xest,3));
                daoed_est.exist_ind_all(seq_tracking.Ts == t) = ~isnan(squeeze(state_est_mcmc.Xest(1,t,seq_tracking.Ts == t)))';
                daoed_est.notclutter_ind = zeros(size(daoed_est.exist_ind_all)); daoed_est.total_no_of_objects = length(daoed_est.exist_ind_all); daoed_est.misseddetect_counter = zeros(size(daoed_est.exist_ind_all)); daoed_est.detect_counter = zeros(size(daoed_est.exist_ind_all));
                [daoed_est, object_class_est, state_est] = initialize_tracking_usingMCMCest(object_class_est, state_est, daoed_est, data_aug.measurements_all{t}, data_aug, hyper_emp, script_params,0,t);

            else

                state_est.Xest(:,seq_tracking.Te > t) = cat(1, squeeze(state_est_mcmc.Xest(1:2,t,seq_tracking.Te > t)), squeeze(state_est_mcmc.Xest(3:4,t,seq_tracking.Te > t))./1e3);
                state_est.Xest(:,seq_tracking.Te == t) = cat(1, squeeze(state_est_mcmc.Xest(1:2,t,seq_tracking.Te == t)), squeeze(state_est_mcmc.Xest(3:4,t,seq_tracking.Te == t)));
                daoed_est.exist_ind_all(seq_tracking.Ts == t) = ~isnan(squeeze(state_est_mcmc.Xest(1,t,seq_tracking.Ts == t)))'; %daoed_est.detect_counter(seq_tracking.Ts == t) = ~isnan(squeeze(state_est_mcmc.Xest(1,t,seq_tracking.Ts == t)))';

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% data association and object existence detection (DAOED) stage
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                [daoed_est, object_class_est, state_est] = detect_objects_and_associate_data(object_class_est, state_est, daoed_est, data_aug.measurements_all{t}, data_aug, hyper_emp, script_params,0,t);

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% EKF-based state estimation stage
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            state_est = track_objects_usingMCMCest(daoed_est, object_class_est, state_est, squeeze(state_est_mcmc.Xest(:,t,seq_tracking.Te > t)), seq_tracking.Te > t, data_aug.measurements_all{t}, data_aug, hyper_emp, script_params);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% store estimates
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %exists_and_detected_ind{t} = (daoed_est.exist_ind_all > 0) & (daoed_est.notclutter_ind > 0);
            exists_and_detected_ind = [exists_and_detected_ind zeros(1, daoed_est.total_no_of_objects - length(exists_and_detected_ind))];
            exists_and_detected_ind = ((daoed_est.exist_ind_all > 0) & (daoed_est.notclutter_ind > 0)) + exists_and_detected_ind;
            Xest_all{t} = NaN(4, daoed_est.total_no_of_objects);
            Xest_all{t}(:, daoed_est.exist_ind_all > 0) = state_est.Xest(:, daoed_est.exist_ind_all > 0);
            dest_all{t} = NaN(2, daoed_est.total_no_of_objects);
            dest_all{t}(1, daoed_est.exist_ind_all > 0) = object_class_est.Dest(1,1,object_class_est.Cest(daoed_est.exist_ind_all > 0));
            dest_all{t}(2, daoed_est.exist_ind_all > 0) = object_class_est.Dest(2,2,object_class_est.Cest(daoed_est.exist_ind_all > 0));
            qest_all{t} = NaN(2, daoed_est.total_no_of_objects);
            qest_all{t}(1, daoed_est.exist_ind_all > 0) = object_class_est.Qest(1,1,object_class_est.Cest(daoed_est.exist_ind_all > 0));
            qest_all{t}(2, daoed_est.exist_ind_all > 0) = object_class_est.Qest(3,3,object_class_est.Cest(daoed_est.exist_ind_all > 0));
            for n = find(daoed_est.exist_ind_all)
                Z(:,1:sum(daoed_est.DA.' == n),t,n) = data_aug.measurements_all{t}(:, daoed_est.DA.' == n);
            end
            Z(:,1:sum(daoed_est.DA.' == 0),t,end) = data_aug.measurements_all{t}(:, daoed_est.DA.' == 0);
        end
        exists_and_detected_ind = exists_and_detected_ind > 0;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% process estimates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~, no_of_objects_cell] = cellfun(@size, Xest_all);
        max_no_of_objects = max(no_of_objects_cell);
        X = NaN(4, script_params.maxT, max_no_of_objects);
        d = NaN(2, script_params.maxT, max_no_of_objects);
        q = NaN(2, script_params.maxT, max_no_of_objects);
        D = NaN(2,2, script_params.maxT, max_no_of_objects);
        for t = 1:script_params.maxT
            X(:,t,exists_and_detected_ind(1:size(Xest_all{t},2))) = Xest_all{t}(:,exists_and_detected_ind(1:size(Xest_all{t},2)));
            d(:,t,exists_and_detected_ind(1:size(Xest_all{t},2))) = dest_all{t}(:,exists_and_detected_ind(1:size(Xest_all{t},2)));
            q(:,t,exists_and_detected_ind(1:size(Xest_all{t},2))) = qest_all{t}(:,exists_and_detected_ind(1:size(Xest_all{t},2)));
            for i = 1:size(Xest_all{t},2)
                D(:,:,t,i) = diag([dest_all{t}(:,i)]);
            end
        end

        %[~, maxNactive] = find(~isnan(squeeze(X(1,:,:))), 1, 'last');
        maxNactive = sum(exists_and_detected_ind);
        nactive = zeros(maxNactive, script_params.maxT);

        Ts = NaN(1, maxNactive);
        Te = NaN(1, maxNactive);
        iprime = 0;
        for i = find(exists_and_detected_ind)
            iprime = iprime + 1;
            nactive(iprime, :) = ~isnan(X(1, :, i));
            if any(find(nactive(iprime, :), 1, 'first'))
                Ts(1, iprime) = find(nactive(iprime, :), 1, 'first');
                Te(1, iprime) = find(nactive(iprime, :), 1, 'last');
            end
        end
        maxTactive = find(sum(nactive, 1), 1, 'last');
        %maxNactive = max_no_of_objects;

        M = NaN(maxTactive, maxNactive);
        nprime = 0;
        for n = find(exists_and_detected_ind)
            nprime = nprime + 1;
            for t = 1:maxTactive
                if all(isnan(Z(1, :, t, n)))
                    M(t, nprime) = 0;
                else
                    M(t, nprime) = find(~isnan(Z(1, :, t, n)), 1, 'last');
                end
            end
        end

        Z = cat(4, Z(:,:,1:script_params.maxT,exists_and_detected_ind), Z(:,:,1:script_params.maxT,end));
        
        seq_tracking = struct('X', X, 'd', d, 'q', q, 'D', D, 'Z', Z, 'object_class_est', object_class_est, 'state_est', state_est, 'daoed_est', daoed_est, 'Ts', Ts, 'Te', Te, 'maxTactive', maxTactive, 'maxNactive', maxNactive, 'nactive', nactive, 'M', M);
    end
end

%results = struct('X', state_est_mcmc.Xest, 'D', object_class_est_mcmc.Dest, 'Z', Z, 'object_class_est', object_class_est_mcmc, 'state_est', state_est_mcmc, 'daoed_est', daoed_est, 'Ts', Ts, 'Te', Te, 'maxTactive', maxTactive, 'maxNactive', maxNactive, 'nactive', nactive, 'M', M);

results = struct('X', state_est_mcmc.Xest, 'D', object_class_est_mcmc.Dest(:,:,object_class_est_mcmc.Cest), 'q', q, 'lambdaM', object_class_est_mcmc.lambdaMest(object_class_est_mcmc.Cest), 'Z', seq_tracking.Z, 'object_class_est', object_class_est_mcmc, 'state_est', state_est_mcmc, 'daoed_est', seq_tracking.daoed_est, 'Ts', seq_tracking.Ts, 'Te', seq_tracking.Te, 'maxTactive', seq_tracking.maxTactive, 'maxNactive', seq_tracking.maxNactive, 'nactive', seq_tracking.nactive, 'M', seq_tracking.M);
end