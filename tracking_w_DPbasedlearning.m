function results = tracking_w_DPbasedlearning(data, hyperparams, script_params, mcmc_params)

exists_and_detected_ind = [];
for t = 1:script_params.maxT

        if t==1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% initialization
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            daoed_est = initialize_tracking(data, hyperparams, script_params);
            state_est = struct('Xest', daoed_est.X, 'Sigmax', repmat(script_params.SigmaXinit, 1, 1, daoed_est.no_of_objectst), 'rotangle', zeros(1,daoed_est.no_of_objectst), 'Sigmarot', repmat(script_params.Sigmarot0, 1, 1, daoed_est.no_of_objectst), 'driving_noise_est', zeros(4, daoed_est.no_of_objectst) ,'no_of_objectstminus1', 0);
            %state_est = struct('Xest', zeros(4, daoed_est.no_of_objectst), 'Sigmax', repmat(SigmaXinit, 1, 1, daoed_est.no_of_objectst), 'rotangle', zeros(1,daoed_est.no_of_objectst), 'Sigmarot', repmat(Sigmarot0, 1, 1, daoed_est.no_of_objectst), 'driving_noise_est', zeros(4, daoed_est.no_of_objectst) ,'no_of_objectstminus1', 0);
            object_class_est = struct('Dest', repmat(hyperparams.D0, 1, 1, daoed_est.no_of_objectst), 'Qest', repmat(hyperparams.Sigmav0, 1, 1, daoed_est.no_of_objectst), 'Cest',(1:daoed_est.no_of_objectst).', 'existsteps', zeros(1, daoed_est.no_of_objectst), 'Msum', zeros(1, daoed_est.no_of_objectst), 'Zsqsum', zeros(2, daoed_est.no_of_objectst), 'Vsqsum', zeros(4, daoed_est.no_of_objectst));


        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% data association and object existence detection (DAOED) stage
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [daoed_est, object_class_est, state_est] = detect_objects_and_associate_data_FUSION(object_class_est, state_est, daoed_est, data.measurements_all{t}, data, hyperparams, script_params,0,t);

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% EKF-based state estimation stage
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        state_est = track_objects(daoed_est, object_class_est, state_est, data.measurements_all{t}, data, hyperparams, script_params);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% DP-based clustering stage
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        object_class_est = learn_objectclasses(state_est, daoed_est, object_class_est, data.measurements_all{t}, data, hyperparams, script_params, mcmc_params, t);

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
            R = calcrotmat(X(:, t, i));
            D(:,:,t,i) = R * diag([dest_all{t}(:,i)]) * R.';
        end
    end

    nactive = zeros(max_no_of_objects, script_params.maxT);
    % [~, maxNactive] = find(~isnan(squeeze(Xdpm(1,:,:))), 1, 'last');

    Ts = NaN(1, max_no_of_objects);
    Te = NaN(1, max_no_of_objects);
    for i = 1:max_no_of_objects
        nactive(i, :) = ~isnan(X(1, :, i));
        if any(find(nactive(i, :), 1, 'first'))
            Ts(1, i) = find(nactive(i, :), 1, 'first');
            Te(1, i) = find(nactive(i, :), 1, 'last');
        end
    end
    maxTactive = find(sum(nactive, 1), 1, 'last');
    maxNactive = max_no_of_objects;


results = struct('X', X, 'd', d, 'q', q, 'D', D, 'object_class_est', object_class_est, 'state_est', state_est, 'daoed_est', daoed_est, 'Ts', Ts, 'Te', Te, 'maxTactive', maxTactive, 'maxNactive', maxNactive, 'nactive', nactive);
end