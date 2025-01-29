function [initialization, mcmc_ind] = initialize_tracking_alt(data, hyper, script, X_init, Drot, lambdaMest)

%use symmetric extents based on D0 average
%davg = 1/2*(hyper.D0(1,1) + hyper.D0(2,2));
%Davg = [davg, 0; 0, davg];
invDrot = NaN(size(Drot));
for i = 1:size(Drot,3)
    invDrot(:,:,i) = Drot(:,:,i)^-1;
end
%% compute measurement partitions according to em algorithm
[clustered_meas, meas_to_cluster_labels, no_of_newobjects_em, clutterindex] = cluster_measurements_em(data.measurements_all{1}, hyper, script, size(X_init,2), script.cov_comp.*Drot, X_init(1:2, :), lambdaMest);

%% assign measurement partitions to existing object states according to gated NN
cluster_to_object_labels = perform_gated_nearest_neighbor_assignment(clustered_meas(:, 1:clutterindex), X_init, Drot, script.thresh);

%% assign measurements to objects according to assigned measurement partitions
meas_to_object_labels = mod(meas_to_cluster_labels, clutterindex);
meas_to_object_labels = zeros(size(meas_to_cluster_labels));
existing_object_ind = 1:size(X_init,2);
for k = find(cluster_to_object_labels > 0)
    %meas_to_object_labels(meas_to_cluster_labels == k) = find(cluster_to_object_labels == k, 1, 'first');
    meas_to_object_labels(meas_to_cluster_labels == k) = existing_object_ind(cluster_to_object_labels(k));
end
if any(meas_to_object_labels > 0)
    [meas_to_object_labels(meas_to_object_labels > 0)] = perform_daod_gating(X_init(1:2,:), meas_to_object_labels(meas_to_object_labels > 0), data.measurements_all{1}(:, meas_to_object_labels > 0), script.initial_thresh, invDrot);
end

%em algorithm for computing initial DA, number of objects, and states
%[X_em_loc, ~, DA1] = initial_state_estimate_em(data.measurements_all{1}, hyper, script);
%[X_em_loc, ~, DA1] = initial_em_known_no_objects(data.measurements_all{1}, hyper, script, sum(data.nactive(:,1)));

%gating procedure to assign measurements to clutter
%[DA1] = perform_init_gating(X_em_loc, DA1, data.measurements_all{1}, script.initial_thresh, invDavg);

%calculate estimate of initial state
no_of_objects_update = sum(unique(meas_to_object_labels) > 0);
%X_init_updated = NaN(4, no_of_objects_update);
[uniq_object_ind, ~, ~] = unique(meas_to_object_labels(meas_to_object_labels > 0));
X_init_updated = X_init(:, uniq_object_ind);
mcmc_ind = zeros(1, size(X_init,2));
mcmc_ind(uniq_object_ind) = find(uniq_object_ind);

%calculate existence indicator and update DA
exist_ind = ~isnan(X_init_updated(1,:));
notclutter_ind = exist_ind;
[~, ~, meas_to_object_labels(meas_to_object_labels > 0)] = unique(meas_to_object_labels(meas_to_object_labels > 0));
meas_to_object_labels = meas_to_object_labels;

initialization = struct("DA", meas_to_object_labels, "exist_indt", exist_ind, "notclutter_ind", notclutter_ind, "X", X_init_updated, "no_of_objectst", no_of_objects_update, "total_no_of_objects", no_of_objects_update, 'misseddetect_counter', exist_ind .* 0, 'detect_counter', exist_ind .* 0, 'exist_ind_all', exist_ind, 'no_of_legacyobjectst', 0, 'no_of_newobjectst', no_of_objects_update);
end