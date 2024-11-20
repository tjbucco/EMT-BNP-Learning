function [DA_updated] = perform_daod_gating(X_loc, DA, measurements, gating_threshold, invCov)

nu_vec = measurements - X_loc(1:2, DA);
nu_vec = nu_vec.';
mahalanobis_dist = NaN(size(nu_vec, 1), 1);
for m = 1:length(DA)
    mahalanobis_dist(m) = (transpose(nu_vec(m,:).')*invCov(:,:,DA(m))*nu_vec(m,:).')^2;
end
%mahalanobis_dist = sum((nu_vec*invCov.').*nu_vec, 2).^2;%mahalanobis_dist = (transpose(nu)*invCov*nu)^2;
outside_gate = mahalanobis_dist > gating_threshold;
outside_gate(isnan(mahalanobis_dist)) = 1;
DA_updated = DA;
DA_updated(outside_gate) = 0;

end