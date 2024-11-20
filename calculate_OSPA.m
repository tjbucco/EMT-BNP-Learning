function [OSPA_metric, OSPA_per_time, loc_error_per_time, card_error_per_time, mse_parameters, mse_parameters_per_time, match_percentage_all] = calculate_OSPA(X, Xest, Ts, Te, Ts_est, Te_est, p, cutoff, alpha, cutoff_mse, D, Dest, q, qest, lambda, lambdaest, C, Cest)

%% assign estimated tracks to labels corresponding with true tracks
assignment_matrix = associate_esttracks_w_truetracks(X, Xest); %omega_tracks, column index corresponds to true track, row to estimated track
OSPA_metric = 0;
mse_parameters = zeros(5,1);
OSPA_per_time = zeros(1, size(X,2));
mse_parameters_per_time = zeros(5,size(X,2));
match_percentage_all = 0;
active_total = 0;

% if nargin > 15
%     [opt_map_Ct] = match_classes(C, Cest, D, Dest, q, qest);
% end

for t = 1:size(X,2)


    [X_cardinality] = sum(~isnan(squeeze(X(1,t,:))).');
    [Xest_cardinality] = sum(~isnan(squeeze(Xest(1,t,:))).');
    ntactive = t >= Ts & t <= Te;
    ntactive_est = t >= Ts_est & t <= Te_est;
    active_total = any(ntactive_est) + active_total;
    [X_cardinality] = size(find(ntactive),2);
    [Xest_cardinality] = size(find(ntactive_est),2);
    X_card_larger = max(X_cardinality, Xest_cardinality);

    match_percentage= 0;
    error_parameters = zeros(5, 1);
    if X_cardinality == 0 && Xest_cardinality == 0
        error_location = 0;
        error_cardinality = 0;
        error_tracklabel = 0;
        error_parameters = zeros(5, 1);
    else
        %% Calculate localization error, i.e., the error in the target states' position/location
        if ~isempty(X_cardinality) && ~isempty(Xest_cardinality)
            if X_card_larger == X_cardinality
                if nargin > 10
                    %[error_location, opt_map] = calculate_GWD_error(squeeze(X(:, t, :)), squeeze(Xest(:, t, :)), X_cardinality, Xest_cardinality, X_card_larger, ntactive, ntactive_est, cutoff, p, D, Dest);
                    [error_location, opt_map] = calculate_GWD_error(squeeze(X(:, t, :)), squeeze(Xest(:, t, :)), X_cardinality, Xest_cardinality, X_card_larger, ntactive, ntactive_est, cutoff, p, D, Dest);
                    [error_parameters] = calculate_parameter_MSE(squeeze(X(:, t, :)), squeeze(Xest(:, t, :)), X_cardinality, Xest_cardinality, X_card_larger, ntactive, ntactive_est, cutoff_mse, p, D, Dest, q, qest, lambda, lambdaest);
                    %[error_parameters] = calculate_parameter_MSE(squeeze(X(:, t, :)), squeeze(Xest(:, t, :)), X_cardinality, Xest_cardinality, X_card_larger, ntactive, ntactive_est, cutoff, p, D, Dest, q, qest);
                    if nargin > 15
                        [match_percentage] = calculate_classification_matches(squeeze(X(:, t, :)), squeeze(Xest(:, t, :)), X_cardinality, Xest_cardinality, X_card_larger, ntactive, ntactive_est, C, Cest);
                    end
                else
                    [error_location, opt_map] = calculate_local_error(squeeze(X(:, t, :)), squeeze(Xest(:, t, :)), X_cardinality, Xest_cardinality, X_card_larger, ntactive, ntactive_est, cutoff, p);
                end
            else
                if nargin > 10
                    %[error_location, opt_map] = calculate_GWD_error(squeeze(Xest(:, t, :)), squeeze(X(:, t, :)), Xest_cardinality, X_cardinality, X_card_larger, ntactive_est, ntactive, cutoff, p, Dest, D);
                    [error_location, opt_map] = calculate_GWD_error(squeeze(Xest(:, t, :)), squeeze(X(:, t, :)), Xest_cardinality, X_cardinality, X_card_larger, ntactive_est, ntactive, cutoff, p, Dest, D);
                    [error_parameters] = calculate_parameter_MSE(squeeze(Xest(:, t, :)), squeeze(X(:, t, :)), Xest_cardinality, X_cardinality, X_card_larger, ntactive_est, ntactive, cutoff_mse, p, Dest, D, qest, q, lambdaest, lambda);
                    %[error_parameters] = calculate_parameter_MSE(squeeze(Xest(:, t, :)), squeeze(X(:, t, :)), Xest_cardinality, X_cardinality, X_card_larger, ntactive_est, ntactive, cutoff, p, Dest, D, qest, q);
                    if nargin > 15
                        [match_percentage] = calculate_classification_matches(squeeze(Xest(:, t, :)), squeeze(X(:, t, :)), Xest_cardinality, X_cardinality, X_card_larger, ntactive_est, ntactive, Cest, C);
                    end
                else
                    [error_location, opt_map] = calculate_local_error(squeeze(Xest(:, t, :)), squeeze(X(:, t, :)), Xest_cardinality, X_cardinality, X_card_larger, ntactive_est, ntactive, cutoff, p);
                end
            end
        else
            error_location = 0;
        end
        %% Calculate CARDINALITY error, i.e., the incorrect number of inferred targets over time
        error_cardinality = abs(X_cardinality - Xest_cardinality) * cutoff^p / X_card_larger;

        %% Calculate LABELING error
        if X_card_larger == X_cardinality
            error_tracklabel = calculate_tracklabeling_error(X, Xest, X_cardinality, Xest_cardinality, ntactive, ntactive_est, opt_map, assignment_matrix);
            %error_dataassociation = calculate_dataassociation_errors(da_true(:,t), da_est_remapped(:,t));
        else
            error_tracklabel = calculate_tracklabeling_error(Xest, X, Xest_cardinality, X_cardinality, ntactive_est, ntactive, opt_map, assignment_matrix.');
            %error_dataassociation = calculate_dataassociation_errors(da_est_remapped(:,t), da_true(:,t));
        end
        error_tracklabel = alpha ^ p / X_card_larger * error_tracklabel;
        %error_dataassociation = alpha ^ p / X_card_larger * error_dataassociation;
    end
    %% Calculate OSPA (Optimal Subpattern Metric)
    loc_error_per_time(t) = error_location;
    card_error_per_time(t) = error_cardinality;
    OSPA_per_time(t) = (error_location + error_cardinality + error_tracklabel) ^ (1 / p);
    OSPA_metric = (error_location + error_cardinality + error_tracklabel) ^ (1 / p) + OSPA_metric;
    %% Calculate MSE (mean squared error)
    mse_parameters_per_time(:,t) = (error_parameters + mse_parameters)./(active_total);
    mse_parameters = error_parameters + mse_parameters;
    %% Calculate classification results
    match_percentage_all = match_percentage_all + match_percentage;
end
OSPA_metric = OSPA_metric/(active_total);
mse_parameters = mse_parameters./(active_total);
match_percentage_all = match_percentage_all/(active_total);
end