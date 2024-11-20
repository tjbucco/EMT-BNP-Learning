function [total_param_e, opt_map] = calculate_parameter_MSE(X1, X2, X1_card, X2_card, X_card_larger, n1active, n2active, cutoff, p, D1, D2, q1, q2, l1, l2)

error_matrix = NaN(X1_card, X2_card);
X1_labels = find(n1active);
X2_labels = find(n2active);
for i = 1:X1_card
    for j = 1:X2_card
        jointly_exist = ~isnan(X1(1,X1_labels(i))) & ~isnan(X2(1,X2_labels(j)));
        if any(jointly_exist)
            euclidean_distance = norm(X1(1:2, X1_labels(i)) - X2(1:2, X2_labels(j)));
            extent_serror = diag(D1(:, :, X1_labels(i)).^2 + D2(:, :, X2_labels(j)).^2 - 2 * D1(:, :, X1_labels(i)) .* D2(:, :, X2_labels(j)));
            driving_noise_var_serror =  q1(:, X1_labels(i)).^2 + q2(:, X2_labels(j)).^2 - 2 * q1(:, X1_labels(i)) .* q2(:, X2_labels(j));
            measurement_rate_serror = l1(X1_labels(i)).^2 + l2(X2_labels(j)).^2 - 2 * l1(X1_labels(i)) .* l2(X2_labels(j));
        else
            euclidean_distance = Inf;
            extent_serror = Inf;
            driving_noise_var_serror = Inf;
            measurement_rate_serror = Inf;
        end
        error_matrix(i,j) = euclidean_distance;
    end
end
[error_location, min_idx] = min(error_matrix(:,:), [], 1);
opt_map = X1_labels(min_idx);
extent_serror = [0;0]; driving_noise_var_serror = [0;0]; measurement_rate_serror = 0;
if X2_card > 0
    for i = 1:X2_card
        extent_serror = extent_serror + min(cutoff(1:2), diag(D1(:, :, opt_map(i)).^2 + D2(:, :, X2_labels(i)).^2 - 2 * D1(:, :, opt_map(i)) .* D2(:, :, X2_labels(i))));
        driving_noise_var_serror =  driving_noise_var_serror +  min(cutoff(3:4), q1(:, opt_map(i)).^2 + q2(:, X2_labels(i)).^2 - 2 * q1(:, opt_map(i)) .* q2(:, X2_labels(i)));
        measurement_rate_serror = measurement_rate_serror + min(cutoff(5), l1(opt_map(i)).^2 + l2(X2_labels(i)).^2 - 2 * l1(opt_map(i)) .* l2(X2_labels(i)));
    end
total_param_e = [extent_serror;driving_noise_var_serror;measurement_rate_serror]./ (1*X2_card);
else
    total_param_e = [extent_serror;driving_noise_var_serror;measurement_rate_serror];
end
%total_param_e = min(cutoff, [extent_serror;driving_noise_var_serror])./ X2_card;
end