function [daoed_est, object_class_est_update, state_est_update] = detect_objects_and_associate_data(object_class_est, state_est, daoed_est_prev, measurements, data, hyper, script, knownparams_flag, t)

if nargin < 8
    knownparams_flag = 0;
end

if isempty(measurements)
    exist_ind_all = daoed_est_prev.exist_ind_all;
    notclutter_ind = daoed_est_prev.notclutter_ind;
    misseddetect_counter = zeros(1, length(exist_ind_all));
    detect_counter = zeros(1, length(exist_ind_all));
    for i = find(daoed_est_prev.exist_ind_all)
        misseddetect_counter(i) = daoed_est_prev.misseddetect_counter(i) + 1;
        if misseddetect_counter(i) >= script.death_heuristic_memory
            exist_ind_all(i) = 0;
        end
    end
    total_no_of_objects = length(exist_ind_all);
    no_of_legacyobjects = sum(exist_ind_all(1:daoed_est_prev.total_no_of_objects));
    daoed_est = daoed_est_prev;
    daoed_est.DA = []; daoed_est.no_of_legacyobjects = no_of_legacyobjects;
    daoed_est.total_no_of_objects = total_no_of_objects; daoed_est.misseddetect_counter = misseddetect_counter;
    daoed_est.exist_ind_all = exist_ind_all; daoed_est.notclutter_ind = notclutter_ind; daoed_est.detect_counter = detect_counter;
end

object_class_est_update = object_class_est;
state_est_update = state_est;
X_new = [];
X_predict = NaN(size(state_est.Xest));
%% compute prediction of object states at next time step
exist_notclutter_ind = (daoed_est_prev.exist_ind_all > 0) & (daoed_est_prev.notclutter_ind > 0);
X_predict(:, daoed_est_prev.exist_ind_all > 0) = hyper.F * state_est.Xest(:, daoed_est_prev.exist_ind_all > 0);
no_of_legacyobjects = sum(daoed_est_prev.exist_ind_all);

if no_of_legacyobjects == 0
    %use symmetric extents based on D0 average
    davg = 1/2*(hyper.D0(1,1) + hyper.D0(2,2));
    Drot = [davg, 0; 0, davg];
    invDavg = Drot^-1;

    [X_em_loc, ~, DA] = initial_state_estimate_em(measurements, hyper, script);
    %[X_em_loc, ~, DA] = initial_em_known_no_objects(measurements, hyper, script, sum(data.nactive(:,t)));

    %gating procedure to assign measurements to clutter
    [DA(DA > 0)] = perform_init_gating(X_em_loc, DA(DA > 0), measurements(DA > 0), script.initial_thresh, invDavg);
    
    %calculate estimate of initial state
    no_of_objects_update = sum(unique(DA) > 0);
    X_new = NaN(4, no_of_objects_update);
    for k = 1:no_of_objects_update
        kprime = unique(DA(DA > 0));
        X_new(1:2,k) = mean(measurements(:,DA == kprime(k)),2) + rand(2, 1,'double');
        %X_em1(1:2,k) = mean(data.measurements_all{1}(:,DA == kprime(k)),2) + em_random_scale.*rand(2, 1 ,'double');
    end
    X_new(1:2, :) = X_new(1:2, :) + rand(2, no_of_objects_update ,'double');
    X_new(3:4, :) = 0;

    %calculate existence indicator and update DA
    newobject_exists = ones(1, no_of_objects_update);
    [~, ~, DA(DA > 0)] = unique(DA(DA > 0));
    meas_to_object_labels = DA + daoed_est_prev.total_no_of_objects;


    if knownparams_flag
        object_class_est_update.Cest = [object_class_est_update.Cest.' ones(sum(newobject_exists), 1).'].';
    else
        [D, Q, lambda, C_new] = assignclassCRPheuristic(no_of_objects_update, object_class_est.Dest(:,:,~isnan(object_class_est.Dest(1,1,:))), object_class_est.Qest(:,:,~isnan(object_class_est.Qest(1,1,:))), object_class_est.lambdaMest(~isnan(object_class_est.lambdaMest)), object_class_est.Cest, hyper);
        if isempty(object_class_est.Cest)
            object_class_est_update.Cest = C_new(newobject_exists > 0);
            object_class_est_update.Dest = D(:,:, unique(C_new(newobject_exists > 0)));
            object_class_est_update.Qest = Q(:,:, unique(C_new(newobject_exists > 0)));
            object_class_est_update.lambdaMest = lambda(unique(C_new(newobject_exists > 0)));
        else
            object_class_est_update.Cest = [object_class_est_update.Cest; C_new(newobject_exists > 0)];
            object_class_est_update.Dest = cat(3, object_class_est_update.Dest(:,:,~isnan(object_class_est_update.Dest(1,1,:))), D(:,:, unique(C_new(((C_new > max(object_class_est.Cest)).' .* (newobject_exists)) > 0))));
            object_class_est_update.Qest = cat(3, object_class_est_update.Qest(:,:,~isnan(object_class_est_update.Qest(1,1,:))), Q(:,:, unique(C_new(((C_new > max(object_class_est.Cest)).' .* (newobject_exists)) > 0))));
            object_class_est_update.lambdaMest = [object_class_est_update.lambdaMest(~isnan(object_class_est_update.lambdaMest)); lambda(unique(C_new(((C_new > max(object_class_est.Cest)).' .* (newobject_exists)) > 0)))];
        end
    end
    object_class_est_update.existsteps = [object_class_est_update.existsteps zeros(1, no_of_objects_update)];
    object_class_est_update.Msum = [object_class_est_update.Msum zeros(1, no_of_objects_update)];
    object_class_est_update.Zsqsum = [object_class_est_update.Zsqsum zeros(2, no_of_objects_update)];
    object_class_est_update.Vsqsum = [object_class_est_update.Vsqsum zeros(4, no_of_objects_update)];
    state_est_update.Xest = cat(2, state_est.Xest, X_new(:, newobject_exists));
    state_est_update.Sigmax = cat(3, state_est_update.Sigmax, repmat(hyper.Sigmax1, 1, 1, no_of_objects_update));
    state_est_update.rotangle = cat(2, state_est_update.rotangle, zeros(1, no_of_objects_update));
    state_est_update.Sigmarot = cat(3, state_est_update.Sigmarot, repmat(diag([0.05, 0.05, 0.006]), 1, 1, no_of_objects_update));
    state_est_update.driving_noise_est = cat(2, state_est_update.driving_noise_est, zeros(4, no_of_objects_update));
    
    no_of_newobjects_est = no_of_objects_update;

else
    %% rotate Dest for each target according to state_est.a
    object_class_est.Dest(:, :, object_class_est.Cest(exist_notclutter_ind));
    Drot = NaN(2,2,daoed_est_prev.total_no_of_objects);
    invDavg = NaN(2,2,daoed_est_prev.total_no_of_objects);
    nprime = 0;
    for n = find(~isnan(X_predict(1, :)))
        nprime = nprime + 1;
        davg = 1/2*(object_class_est.Dest(1, 1, object_class_est.Cest(n)) + object_class_est.Dest(2, 2, object_class_est.Cest(n)));
        R = calcrotmat(X_predict(:, n));
        Drot(:,:,n) = R * object_class_est.Dest(:, :, object_class_est.Cest(n)) * R.';
        Drot(:,:,n) = davg * eye(2);
        invDavg(:,:,n) = (Drot(:,:,n))^-1;
    end

    %% compute measurement partitions according to em algorithm
    [clustered_meas, meas_to_cluster_labels, no_of_newobjects_em, clutterindex] = cluster_measurements_em(measurements, hyper, script, sum(daoed_est_prev.exist_ind_all > 0), script.cov_comp.*Drot(:,:,daoed_est_prev.exist_ind_all > 0), X_predict(1:2, daoed_est_prev.exist_ind_all > 0), object_class_est.lambdaMest(object_class_est.Cest(daoed_est_prev.exist_ind_all > 0)));

    %% assign measurement partitions to existing object states according to gated NN
    cluster_to_object_labels = perform_gated_nearest_neighbor_assignment(clustered_meas(:, 1:clutterindex), X_predict(:,daoed_est_prev.exist_ind_all > 0), object_class_est.Dest(:, :, object_class_est.Cest(daoed_est_prev.exist_ind_all > 0)), script.thresh);

    %% assign measurements to objects according to assigned measurement partitions
    meas_to_object_labels = mod(meas_to_cluster_labels, clutterindex);
    meas_to_object_labels = zeros(size(meas_to_cluster_labels));
    %for k = unique(cluster_to_object_labels)
    existing_object_ind = find(daoed_est_prev.exist_ind_all);
    for k = find(cluster_to_object_labels > 0)
        %meas_to_object_labels(meas_to_cluster_labels == k) = find(cluster_to_object_labels == k, 1, 'first');
        meas_to_object_labels(meas_to_cluster_labels == k) = existing_object_ind(cluster_to_object_labels(k));
    end
    if any(meas_to_object_labels > 0)
        [meas_to_object_labels(meas_to_object_labels > 0)] = perform_daod_gating(X_predict(1:2,:), meas_to_object_labels(meas_to_object_labels > 0), measurements(:, meas_to_object_labels > 0), script.thresh, invDavg);
    end
    %% birth heuristic
    idx=ismember(cluster_to_object_labels, 1:clutterindex); % idx
    unassigned_clusters = clustered_meas(:, ~idx);
    unassigned_clusters = unassigned_clusters(:, ~isnan(unassigned_clusters(1,:)));
    no_of_newobjects_est = min(sum(~isnan(unassigned_clusters(1,:))), script.Nbirth);
    if no_of_newobjects_est > 0
        if knownparams_flag
            C_new = ones(no_of_newobjects_est,1);
            D = object_class_est.Dest;
            lambda = object_class_est.lambdaMest;
        else
            [D, Q, lambda, C_new] = assignclassCRPheuristic(no_of_newobjects_est, object_class_est.Dest(:,:,~isnan(object_class_est.Dest(1,1,:))), object_class_est.Qest(:,:,~isnan(object_class_est.Qest(1,1,:))), object_class_est.lambdaMest(~isnan(object_class_est.lambdaMest)), object_class_est.Cest, hyper);
        end
        invDnew = NaN(2,2,no_of_newobjects_est);
        for n = 1:no_of_newobjects_est
            davg = 1/2*(D(1, 1, C_new(n)) + D(2, 2, C_new(n)));
            invDnew(:,:,n) = (davg * eye(2))^-1;
        end
        X_new = NaN(4, no_of_newobjects_est);
        %[X_new_loc, meas_to_newbirth_labels, no_of_newobjects_em, clutterindex] = cluster_measurements_em(measurements(meas_to_object_labels == 0), hyper, script, 1, Davg);
        %[X_new_loc, ~, DA1] = initial_em_known_no_objects(measurements(meas_to_object_labels == 0), hyper, script, no_of_newobjects_est);
        %X_new(1:2, :) = clustered_meas(:, object_to_cluster_labels == 0);
        X_new(1:2,:) = unassigned_clusters(1:2,1:no_of_newobjects_est) + (rand(1, 2 ,'double')).';
        X_new(3:4, :) = 0;

        cluster_to_objects_and_newobjects_labels = perform_gated_nearest_neighbor_assignment(clustered_meas(:, 1:clutterindex), cat(2, X_predict(:,daoed_est_prev.exist_ind_all > 0), X_new), cat(3, object_class_est.Dest(:, :, object_class_est.Cest(daoed_est_prev.exist_ind_all > 0)), D(:,:,C_new)), script.thresh_cluster);

        %assign measurements to objects according to assigned measurement partitions
        meas_to_object_labels = mod(meas_to_cluster_labels, clutterindex);
        meas_to_object_labels = zeros(size(meas_to_cluster_labels));
        existing_object_ind = find([daoed_est_prev.exist_ind_all  ones(1, no_of_newobjects_est)]);
        for k = find(cluster_to_objects_and_newobjects_labels > 0)
            meas_to_object_labels(meas_to_cluster_labels == k) = existing_object_ind(cluster_to_objects_and_newobjects_labels(k));
        end
        if any(meas_to_object_labels > 0)
            [meas_to_object_labels(meas_to_object_labels > 0)] = perform_daod_gating([X_predict(1:2,:) X_new(1:2,:)], meas_to_object_labels(meas_to_object_labels > 0), measurements(:, meas_to_object_labels > 0), script.thresh, cat(3, invDavg, invDnew));
        end

        no_of_newobjects_est_round2 = sum(unique(meas_to_object_labels) > daoed_est_prev.total_no_of_objects);
        if no_of_newobjects_est_round2 > 0
            newobject_exists = zeros(1, no_of_newobjects_est);
            iprime = 0;
            for i = unique(meas_to_object_labels(meas_to_object_labels > daoed_est_prev.total_no_of_objects))
                iprime = i - daoed_est_prev.total_no_of_objects;
                newobject_exists(iprime) = any(meas_to_object_labels == i);
            end
            
            newobject_exists = newobject_exists > 0;
            [~, ~, temp] = unique(meas_to_object_labels(meas_to_object_labels > daoed_est_prev.total_no_of_objects));
            meas_to_object_labels(meas_to_object_labels > daoed_est_prev.total_no_of_objects) = temp.' + daoed_est_prev.total_no_of_objects;

            temp = C_new(newobject_exists);
            if any(temp(temp > max(object_class_est_update.Cest)))
                temp2 = unique(temp(temp > max(object_class_est_update.Cest)));
                [~, ~, addCnew] = unique(temp2);
                C_new_new = addCnew.' + max(object_class_est_update.Cest);
                D(:,:,C_new_new) = D(:,:, temp2);
                Q(:,:,C_new_new) = Q(:,:, temp2);
                lambda(C_new_new) = lambda(temp2);
                iprime = 0;
                for i = C_new_new
                    iprime = iprime + 1;
                    C_new(C_new == temp2(iprime)) = i;
                end
                %C_new(ismember(C_new, temp(temp > max(object_class_est_update.Cest)))) = C_new_new;
            end

            if knownparams_flag
                object_class_est_update.Cest = [object_class_est_update.Cest ; ones(sum(newobject_exists), 1)];
            else
                object_class_est_update.Cest = [object_class_est_update.Cest ; C_new(newobject_exists)];
                object_class_est_update.Dest = cat(3, object_class_est_update.Dest(:,:,~isnan(object_class_est_update.Dest(1,1,:))), D(:,:, unique(C_new(((C_new > max(object_class_est.Cest)).' .*(newobject_exists)) > 0))));
                object_class_est_update.Qest = cat(3, object_class_est_update.Qest(:,:,~isnan(object_class_est_update.Qest(1,1,:))), Q(:,:, unique(C_new(((C_new > max(object_class_est.Cest)).' .*(newobject_exists)) > 0))));
                object_class_est_update.lambdaMest = [object_class_est_update.lambdaMest(~isnan(object_class_est_update.lambdaMest)); lambda(unique(C_new(((C_new > max(object_class_est.Cest)).' .*(newobject_exists)) > 0)))];
            end
            object_class_est_update.existsteps = [object_class_est_update.existsteps zeros(1, no_of_newobjects_est_round2)];
            object_class_est_update.Msum = [object_class_est_update.Msum zeros(1, no_of_newobjects_est_round2)];
            object_class_est_update.Zsqsum = [object_class_est_update.Zsqsum zeros(2, no_of_newobjects_est_round2)];
            object_class_est_update.Vsqsum = [object_class_est_update.Vsqsum zeros(4, no_of_newobjects_est_round2)];
            state_est_update.Xest = cat(2, state_est.Xest, X_new(:, newobject_exists));
            state_est_update.Sigmax = cat(3, state_est_update.Sigmax, repmat(hyper.Sigmax1, 1, 1, no_of_newobjects_est_round2));
            state_est_update.rotangle = cat(2, state_est_update.rotangle, zeros(1, no_of_newobjects_est_round2));
            state_est_update.Sigmarot = cat(3, state_est_update.Sigmarot, repmat(diag([0.05, 0.05, 0.006]), 1, 1, no_of_newobjects_est_round2));
            state_est_update.driving_noise_est = cat(2, state_est_update.driving_noise_est, zeros(4, no_of_newobjects_est_round2));

        end
        no_of_newobjects_est = no_of_newobjects_est_round2;

    end
end

exist_indt = ones(1, no_of_legacyobjects + no_of_newobjects_est);
exist_ind_all = [daoed_est_prev.exist_ind_all ones(1, no_of_newobjects_est)];
notclutter_ind = [daoed_est_prev.notclutter_ind zeros(1, no_of_newobjects_est)];
misseddetect_counter = zeros(1, length(exist_ind_all));
detect_counter = zeros(1, length(exist_ind_all));
for i = find(exist_ind_all)
    if ~any(meas_to_object_labels == i)
        misseddetect_counter(i) = daoed_est_prev.misseddetect_counter(i) + 1;
        if misseddetect_counter(i) >= script.death_heuristic_memory
            exist_ind_all(i) = 0;
        end
    else
        if i > daoed_est_prev.total_no_of_objects
            detect_counter(i) = detect_counter(i) + 1;
        else
            detect_counter(i) = daoed_est_prev.detect_counter(i) + 1;
        end
        if detect_counter(i) >= script.birth_heuristic_memory
            notclutter_ind(i) = 1;
        end        
    end
end

no_of_objectst = sum(exist_ind_all);
%misseddetect_counter = [misseddetect_counter zeros(1, no_of_newobjects_est)];

total_no_of_objects = length(exist_ind_all);
no_of_legacyobjects = sum(exist_ind_all(1:daoed_est_prev.total_no_of_objects));
%meas_to_object_labels(meas_to_object_labels > 0) = meas_to_object_labels(meas_to_object_labels > 0) + total_no_of_objects - no_of_objectst;
daoed_est = struct('DA', meas_to_object_labels, 'exist_indt', exist_indt, 'notclutter_ind', notclutter_ind, 'detect_counter', detect_counter, 'X_new', X_new, 'no_of_objectst', no_of_objectst, 'no_of_legacyobjectst', no_of_legacyobjects, 'no_of_newobjectst', no_of_newobjects_est, 'total_no_of_objects', total_no_of_objects, 'misseddetect_counter', misseddetect_counter, 'exist_ind_all', exist_ind_all);
end