function [DA_updated] = perform_init_gating(X_loc, DA, measurements, gating_threshold, invCov)

nu_vec = measurements - X_loc(1:2, DA);
nu_vec = nu_vec.';
mahalanobis_dist = sum((nu_vec*invCov.').*nu_vec, 2).^2;
%mahalanobis_dist = (transpose(nu)*invCov*nu)^2;
outside_gate = mahalanobis_dist > gating_threshold;
outside_gate(isnan(mahalanobis_dist)) = 1;
DA_updated = DA;
DA_updated(outside_gate) = 0;

end