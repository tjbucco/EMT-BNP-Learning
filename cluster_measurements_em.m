function [X_initial, cluster_labels, no_of_new_objects, clutterindex, X_kmeans] = cluster_measurements_em(measurements, hyper, script, no_of_objects, Cov_extent, Xpredict, rate_noofmeas)

no_of_measurements = sum(~isnan(measurements(1,:)),2);
cluster_labels_kmeans = NaN(no_of_measurements, no_of_objects + script.Nbirth);
X_kmeans = NaN(2, no_of_objects + script.Nbirth + 1, script.Nbirth + 1);
X_em = X_kmeans;
X_initial = NaN(2, no_of_objects + script.Nbirth);
dataEM = NaN(3, no_of_measurements, script.Nbirth + 1);
BIC = NaN(1, script.Nbirth + 1);
AIC = NaN(1, script.Nbirth + 1);

davg = 1/2*(hyper.D0(1,1) + hyper.D0(2,2));
Cov_avg = [davg, 0; 0, davg];

if no_of_measurements == 0
    cluster_labels = 0;
else
    for k = 1:script.Nbirth + 1
        if no_of_objects + k == 1
            cluster_labels_kmeans(:, 1) = 1;
            X_kmeans(:, 1, 1) = mean(measurements, 2);
        else
            [cluster_labels_kmeans(:, k), X_kmeans(:, 1:(no_of_objects + k), k)] = kmeans_sorber(measurements, no_of_objects + k);
            %%todo: incorporate predictions into this part. figure out
            %%which X_kmeans match Xpredict and replace those
        end
        X_kmeans(:, 1:(no_of_objects + k - 1), k) = get_nearest_state_munkres(Xpredict, X_kmeans(:, 1:(no_of_objects + k - 1), k));
        points_in_cluster = histc(cluster_labels_kmeans(:, k), 1:no_of_objects + k);
        %points_in_cluster = [points_in_cluster.', 0].';
        % calculate the weights
        %lambda = zeros(no_of_objects + k, 1);
        lambda = [[rate_noofmeas; hyper.ElambdaM * ones(k - 1, 1)]; hyper.lambdaclutter];
        lambda = lambda + points_in_cluster;
        lambda = lambda./sum(lambda);
%         lambda = [hyper.rate_noofmeas * ones(no_of_objects + k - 1, 1); hyper.lambdaclutter];
%         lambda = lambda./sum(lambda);
        Cov_EM = cat(3, Cov_extent, repmat(Cov_avg, 1, 1, k - 1));
        Cov_EM = cat(3, Cov_EM, hyper.clutter_gauss_approx_cov);
        [dataEM(:, :, k), X_em(:, 1:(no_of_objects + k), k), sigma_em, lambda_em] = EM_fixedsig_and_clutter([measurements ; cluster_labels_kmeans(:, k).'], cat(2, X_kmeans(:, 1:(no_of_objects - 1 + k), k), [0;0]), Cov_EM, lambda.');
        est_no_of_objects = sum(~isnan(X_em(1, 1:(no_of_objects + k), k)));
        X_bic = X_em(:,~isnan(X_em(1, 1:(no_of_objects + k), k)),k);
        sigma_bic = sigma_em(:,:,~isnan(X_em(1, 1:(no_of_objects + k), k)));
        loglike = 0;
        likelihood = 0;
        for ii = 1:no_of_measurements
            for kk = 1:est_no_of_objects
                likelihood = lambda_em(1, kk) * mvnpdf_fast(measurements(:,ii), X_bic(:,kk), diag(sigma_bic(:, :, kk)).') + likelihood;
            end
            loglike = log(likelihood) + loglike;
        end
        BIC(1, k) = - loglike + 0.5 * (2 + 4 + 1) * est_no_of_objects * log(no_of_measurements);
        AIC(1, k) = - loglike + 0.5 * (2 + 4 + 1) * est_no_of_objects * 2;
    end
end
[~, model_select] = min(AIC);
maxNactive = sum(~isnan(X_em(1, :, model_select)));
%X_em = X_em(:,~isnan(X_em(1, :, model_select)),model_select);
X_em = X_em(:,:,model_select);
X_initial = X_em;
% X_initial(:, 1:size(X_em, 2)) = X_em;
% X_initial = X_initial(~isnan(X_initial(:,:)));
no_of_new_objects = model_select - 1;
clutterindex = no_of_objects + model_select;
cluster_labels = dataEM(3, :, model_select);
%cluster_labels = cluster_labels(~isnan(cluster_labels));
end