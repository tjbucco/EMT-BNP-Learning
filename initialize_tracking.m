function [initialization] = initialize_tracking(data, hyper, script, em_random_scale)

if nargin < 4
    em_random_scale = 1/1000;
end

%use symmetric extents based on D0 average
davg = 1/2*(hyper.D0(1,1) + hyper.D0(2,2));
Davg = [davg, 0; 0, davg];
invDavg = Davg^-1;


%em algorithm for computing initial DA, number of objects, and states
[X_em_loc, ~, DA1] = initial_state_estimate_em(data.measurements_all{1}, hyper, script, script.use_binoN_init);
% [X_em_loc, ~, DA1] = initial_em_known_no_objects(data.measurements_all{1}, hyper, script, script.binoN);

%gating procedure to assign measurements to clutter
[DA1(DA1 > 0)] = perform_init_gating(X_em_loc, DA1(DA1 > 0), data.measurements_all{1}(DA1 > 0), script.initial_thresh, invDavg);

%calculate estimate of initial state
no_of_objects_update = sum(unique(DA1) > 0);
X_em1 = NaN(4, no_of_objects_update); 
for k = 1:no_of_objects_update
    kprime = unique(DA1(DA1 > 0));
    X_em1(1:2,k) = mean(data.measurements_all{1}(:,DA1 == kprime(k)),2) + em_random_scale.*rand(2, 1 ,'double');
    %X_em1(1:2,k) = mean(data.measurements_all{1}(:,DA1 == kprime(k)),2) + em_random_scale.*rand(2, 1 ,'double');
end
X_em1(1:2, :) = X_em1(1:2, :) + em_random_scale.*rand(2, no_of_objects_update ,'double');
X_em1(3:4, :) = 0 + em_random_scale.*rand(2, no_of_objects_update ,'double');
[~, b, a] = get_nearest_state_munkres(X_em1(1:2,:),squeeze(data.X(1:2,1,:)));
X_em1(3:4, b(b>0)) = data.X(3:4, 1, a);

%calculate existence indicator and update DA
exist_ind = ones(1, no_of_objects_update);
notclutter_ind = exist_ind;
[~, ~, DA1(DA1 > 0)] = unique(DA1(DA1 > 0));
DA1 = DA1;

initialization = struct("DA", DA1, "exist_indt", exist_ind, "notclutter_ind", notclutter_ind, "X", X_em1, "no_of_objectst", no_of_objects_update, "total_no_of_objects", no_of_objects_update, 'misseddetect_counter', exist_ind .* 0, 'detect_counter', exist_ind .* 0, 'exist_ind_all', exist_ind, 'no_of_legacyobjectst', 0, 'no_of_newobjectst', no_of_objects_update);
end