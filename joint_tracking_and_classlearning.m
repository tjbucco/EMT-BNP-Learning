function [results, OSPA_DPclustering, MSE_DPclustering] = joint_tracking_and_classlearning(data, hyperparams, script_params, mcmc_params, seq_tracking, sim)

p = 1; cutoff = 120; cutoff_GWD = cutoff*10^0; cutoff_mse = cutoff*100.*[hyperparams.Elsq_data.^2; hyperparams.Essq_data.^2; hyperparams.ElambdaM_data]; alpha = 2; %OSPA metric parameters

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
    [~, OSPA_DPclustering(:, j)] = calculate_OSPA(data.X, state_est_mcmc.Xest, data.Ts, data.Te, seq_tracking.Ts, seq_tracking.Te, p, cutoff, alpha);
    [~, OSPA_DPclustering_GWD(:, j), MSE_DPclustering(:,j), ~, ~]  = calculate_OSPA(data.X, state_est_mcmc.Xest, data.Ts, data.Te, seq_tracking.Ts, seq_tracking.Te, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), object_class_est_mcmc.Dest(:,:,object_class_est_mcmc.Cest), data.q, q, data.lambdaMstar(data.C), object_class_est_mcmc.lambdaMest(object_class_est_mcmc.Cest), data.C, object_class_est_mcmc.Cest);

    %% insert initial state estimates as measurements
    data_aug = data;
    for i = 1:seq_tracking.maxNactive
        for t = seq_tracking.Ts(i):seq_tracking.Te(i)
            measurements_tent = data_aug.measurements_all{t};
            data_aug.measurements_all{t} = cat(2, measurements_tent, state_est_mcmc.Xest(1:2,t,i));
        end
    end
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
                daoed_est = initialize_tracking(data_aug, hyper_emp, script_params);
                state_est = struct('Xest', daoed_est.X, 'Sigmax', repmat(script_params.SigmaXinit, 1, 1, daoed_est.no_of_objectst), 'rotangle', zeros(1,daoed_est.no_of_objectst), 'Sigmarot', repmat(script_params.Sigmarot0, 1, 1, daoed_est.no_of_objectst), 'driving_noise_est', zeros(4, daoed_est.no_of_objectst) ,'no_of_objectstminus1', 0);
                %state_est = struct('Xest', zeros(4, daoed_est.no_of_objectst), 'Sigmax', repmat(SigmaXinit, 1, 1, daoed_est.no_of_objectst), 'rotangle', zeros(1,daoed_est.no_of_objectst), 'Sigmarot', repmat(Sigmarot0, 1, 1, daoed_est.no_of_objectst), 'driving_noise_est', zeros(4, daoed_est.no_of_objectst) ,'no_of_objectstminus1', 0);
                object_class_est = struct('Dest', repmat(hyper_emp.D0, 1, 1, daoed_est.no_of_objectst), 'Qest', repmat(hyper_emp.Sigmav0, 1, 1, daoed_est.no_of_objectst), 'Cest',zeros(1,daoed_est.no_of_objectst), 'existsteps', zeros(1, daoed_est.no_of_objectst), 'Msum', zeros(1, daoed_est.no_of_objectst), 'Zsqsum', zeros(2, daoed_est.no_of_objectst), 'Vsqsum', zeros(4, daoed_est.no_of_objectst));

                [D, Q, lambda, C_new] = assignclassCRPheuristic(daoed_est.no_of_objectst, object_class_est_mcmc.Dest(:,:,~isnan(object_class_est_mcmc.Dest(1,1,:))), object_class_est_mcmc.Qest(:,:,~isnan(object_class_est_mcmc.Qest(1,1,:))), object_class_est_mcmc.lambdaMest(~isnan(object_class_est_mcmc.lambdaMest)), object_class_est_mcmc.Cest, hyperparams);
                [~, ~, object_class_est.Cest] = unique(C_new);
                object_class_est.Dest = D(:,:, C_new);
                object_class_est.Qest = Q(:,:, C_new);
                object_class_est.lambdaMest = lambda(C_new);

                %object_class_est.existsteps = zeros(1, daoed_est.no_of_objectst);
                %object_class_est.Msum = zeros(1, daoed_est.no_of_objectst);
                %object_class_est.Zsqsum = zeros(2, daoed_est.no_of_objectst);
                %object_class_est.Vsqsum = zeros(4, daoed_est.no_of_objectst);
                
                if any(seq_tracking.Ts == 1)
                    object_class_est = match_objects_to_existingclasses(object_class_est, object_class_est_mcmc, find(seq_tracking.Ts == 1), daoed_est.DA, 1);
                end


            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% data association and object existence detection (DAOED) stage
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [daoed_est, object_class_est, state_est] = detect_objects_and_associate_data(object_class_est, state_est, daoed_est, data_aug.measurements_all{t}, data_aug, hyper_emp, script_params,0,t);

                if any(seq_tracking.Ts == t)
                    object_class_est = match_objects_to_existingclasses(object_class_est, object_class_est_mcmc, find(seq_tracking.Ts == t), daoed_est.DA, t);
                end

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% EKF-based state estimation stage
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            state_est = track_objects(daoed_est, object_class_est, state_est, data_aug.measurements_all{t}, data_aug, hyper_emp, script_params);


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