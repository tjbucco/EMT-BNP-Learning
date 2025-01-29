function [object_class_est] = cluster_objectsEM(state_est, daoed_est, object_class_est_prev, measurements, data, hyper, script, mcmc, t)

%update sufficient statistics
existsteps = object_class_est_prev.existsteps;
Msum = object_class_est_prev.Msum;
Zsqsum = object_class_est_prev.Zsqsum;
Vsqsum = object_class_est_prev.Vsqsum;
 
for i = find(daoed_est.exist_ind_all)
    rot = calcrotmat(state_est.Xest(:,i));
    premeasurements_i = rot.'* (measurements(:, daoed_est.DA == i) - state_est.Xest(1:2,i));
    existsteps(i) = object_class_est_prev.existsteps(i) + 1;
    Msum(i) = object_class_est_prev.Msum(i) + sum(daoed_est.DA == i);
    Zsqsum(:, i) = object_class_est_prev.Zsqsum(:, i) + sum(premeasurements_i.^2, 2);
    Vsqsum(:,i) = object_class_est_prev.Vsqsum(:,i) + sum(state_est.driving_noise_est(:, i).^2, 2);
end
if any(daoed_est.exist_ind_all)
    object_class_est.existsteps = existsteps; object_class_est.Msum = Msum; 
    object_class_est.Zsqsum = Zsqsum; object_class_est.Vsqsum = Vsqsum;
    object_class_est.lambdaMest = object_class_est_prev.lambdaMest;
else
    object_class_est = object_class_est_prev;
    return
end

%%EM algorithm - see He, Valin 2014
no_of_objects = daoed_est.total_no_of_objects;
Cmax = script.Cmax;
cluster_labels_kmeans = NaN(no_of_objects, script.Cmax);
P_kmeans = NaN(6, Cmax, Cmax);
P_em = P_kmeans;
P_initial = NaN(6, Cmax);
dataEM = NaN(7, no_of_objects, Cmax);
BIC = NaN(1, Cmax);
AIC = NaN(1, Cmax);
Dobs = 200*Zsqsum(:, 1:no_of_objects)./Msum(1:no_of_objects);
Dobs(:, Msum(1:no_of_objects) == 0) = 0;
Qobs = 1*Vsqsum(:, 1:no_of_objects)./(existsteps(1:no_of_objects));
%davg = 1/2*(hyper.D0(1,1) + hyper.D0(2,2));
% Cov_avg = [davg, 0; 0, davg];

if no_of_objects == 0
    cluster_labels = 0;
else
    for k = 1:Cmax
        if Cmax == 1
            cluster_labels_kmeans(:, 1) = 1;
            P_kmeans(1:2, 1, 1) = mean(Dobs, 2);
            P_kmeans(3:end, 1, 1) = mean(Qobs, 2);
        else
            [cluster_labels_kmeans(:, k), P_kmeans(:, 1:k, k)] = kmeans_sorber(cat(1, Dobs, Qobs), k);
        end
        points_in_cluster = histc(cluster_labels_kmeans(:, k), unique(cluster_labels_kmeans(:, k)));
        points_in_cluster = [points_in_cluster.', 0].';
        % calculate the weights
        lambda = zeros(k + 1, 1);
        lambda(1:length(points_in_cluster), 1) = points_in_cluster./sum(points_in_cluster);
        Cov_EM = repmat(script.Cov_EM, 1, 1, k);        
        [dataEM(:, :, k), P_em(:, 1:k, k), sigma_em, lambda_em] = EM([cat(1, Dobs, Qobs) ; cluster_labels_kmeans(:, k).'], P_kmeans(:, 1:k, k), Cov_EM, lambda.');
        est_no_of_classes = sum(~isnan(P_em(1, 1:k, k)));
        P_bic = P_em(:,~isnan(P_em(1, 1:k, k)),k);
        sigma_bic = sigma_em(:,:,~isnan(P_em(1, 1:k, k)));
        loglike = 0;
        likelihood = 0;
        for ii = 1:no_of_objects
            for kk = 1:est_no_of_classes
                %likelihood = lambda_em(1, kk) * mvnpdf(log([Dobs(:,ii); Qobs(:,ii)]), log(P_bic(:,kk))) + likelihood;
                likelihood = lambda_em(1, kk) * mvnpdf([Dobs(:,ii); Qobs(:,ii)], P_bic(:,kk)) + likelihood;
            end
            loglike = log(likelihood) + loglike;
        end
        BIC(1, k) = - loglike + 0.5 * (2 + 4 + 1) * est_no_of_classes * log(no_of_objects);
        AIC(1, k) = - loglike + 0.5 * (2 + 4 + 1) * est_no_of_classes * 2;
    end
end
[~, model_select] = min(BIC);
maxNactive = sum(~isnan(P_em(1, :, model_select)));
%P_em = P_em(:,~isnan(P_em(1, :, model_select)),model_select);
P_em = P_em(:,:,model_select);
P_initial = P_em;
% X_initial(:, 1:size(P_em, 2)) = P_em;
% X_initial = X_initial(~isnan(X_initial(:,:)));
cluster_labels = dataEM(end, :, model_select);
%cluster_labels = cluster_labels(~isnan(cluster_labels));
assigned_clusters = unique(cluster_labels);
cprime = 0;
for c = assigned_clusters
    cprime = cprime + 1;
    Dstar(:, :, cprime) = diag(P_em(1:2, c).'.^(1/2));
    Qstar(:, :, cprime) = diag([ones(1,2)*(sum(P_em(3:4, c))/2)^(1/2), ones(1,2)*(sum(P_em(5:6, c))/2)^(1/2)]);
end
[~,~,remap_clusters] = unique(cluster_labels);
object_class_est.Cest = remap_clusters;
object_class_est.Dest = Dstar;
object_class_est.Qest = Qstar;

end