function [X_initial, maxNactive, cluster_labels, no_of_new_objects] = initial_state_estimate_em(measurements, hyper, script, use_binoN)

if nargin < 4
    use_binoN = 0;
end
if use_binoN
    Nbirth = script.binoN;
else
    Nbirth = script.Nbirth;
end

no_of_measurements = sum(~isnan(measurements(1,:)),2);
cluster_labels_kmeans = NaN(no_of_measurements, Nbirth);
X_kmeans = NaN(2, Nbirth + 1, Nbirth + 1);
X_em = X_kmeans;
X_initial = NaN(2, Nbirth);
dataEM = NaN(3, no_of_measurements, Nbirth + 1);
BIC = NaN(1, Nbirth + 1);
AIC = NaN(1, Nbirth + 1);
if Nbirth == 0
    maxNactive = 0;
    cluster_labels = 0;
    return
end

if no_of_measurements == 0
    cluster_labels = 0;
else
    %cluster_labels = kmeans_sorber(measurements,min([no_of_measurements, maxNactive]));
    davg = 1/2*(hyper.D0(1,1) + hyper.D0(2,2));
    Cov_avg = script.cov_comp.*[davg, 0; 0, davg];

    for k = 1:Nbirth + 1
        if k == 1
            cluster_labels_kmeans(:, 1) = 1;
            X_kmeans(:, 1, 1) = mean(measurements, 2);
        else
            [cluster_labels_kmeans(:, k), X_kmeans(:, 1:k, k)] = kmeans_sorber(measurements, k);
        end

        points_in_cluster = histc(cluster_labels_kmeans(:, k), 1:k);
        % calculate the weights
        lambda = [hyper.ElambdaM * ones(k - 1, 1); hyper.lambdaclutter];
        lambda = lambda + points_in_cluster;
        lambda = lambda./sum(lambda);
        Cov_EM = repmat(Cov_avg, 1, 1, k - 1);
        Cov_EM = cat(3, Cov_EM, hyper.clutter_gauss_approx_cov);
        [dataEM(:, :, k), X_em(:, 1:k, k), sigma_em, lambda_em] = EM_fixedsig_and_clutter([measurements ; cluster_labels_kmeans(:, k).'], cat(2,X_kmeans(:, 1:k - 1, k), [0;0]), Cov_EM, lambda.');
        no_of_tar = sum(~isnan(X_em(1, 1:k, k)));
        X_bic = X_em(:,~isnan(X_em(1, 1:k, k)),k);
        sigma_bic = sigma_em(:,:,~isnan(X_em(1, 1:k, k)));
        loglike = 0;
        likelihood = 0;
        for ii = 1:no_of_measurements
            for kk = 1:no_of_tar
%                 likelihood = lambda_em(1, kk) * mvnpdf(measurements(:,ii), X_bic(:,kk), diag(sigma_bic(:, :, kk)).') + likelihood;
                likelihood = lambda_em(1, kk) * mvnpdf_fast(measurements(:,ii), X_bic(:,kk), diag(sigma_bic(:, :, kk)).') + likelihood;
            end
            loglike = log(likelihood) + loglike;
        end
        BIC(1, k) = - loglike + 0.5 * (2 + 4 + 1) * no_of_tar * log(no_of_measurements);
        AIC(1, k) = - loglike + 0.5 * (2 + 4 + 1) * no_of_tar * 2;
    end
end
[~, model_select] = min(AIC);
maxNactive = sum(~isnan(X_em(1, :, model_select)));
%X_em = X_em(:,~isnan(X_em(1, :, model_select)),model_select);
X_em = X_em(:,:,model_select);
%X_initial(:, 1:size(X_em, 2)) = X_em;
X_initial = X_em(:,1:(model_select - 1));
no_of_new_objects = model_select - 1;
clutterindex = model_select;
cluster_labels = dataEM(3, :, model_select);
cluster_labels(cluster_labels==model_select) = 0;

end