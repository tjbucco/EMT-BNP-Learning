function [total_local_e, opt_map] = calculate_GWD_error(X1, X2, X1_card, X2_card, X_card_larger, n1active, n2active, cutoff, p, D1, D2)

error_matrix = NaN(X1_card, X2_card);
X1_labels = find(n1active);
X2_labels = find(n2active);
for i = 1:X1_card
    for j = 1:X2_card
        jointly_exist = ~isnan(X1(1,X1_labels(i))) & ~isnan(X2(1,X2_labels(j)));
        if any(jointly_exist)
            euclidean_distance = norm(X1(1:2, X1_labels(i)) - X2(1:2, X2_labels(j)),2)^2;
            extent_penalty = trace(D1(:, :, X1_labels(i)) + D2(:, :, X2_labels(j)) - 2 * sqrtm(sqrtm(D1(:, :, X1_labels(i))) * D2(:, :, X2_labels(j)) * sqrtm(D1(:, :, X1_labels(i)))));
            gaussian_wasserstein_distance = euclidean_distance + extent_penalty;
        else
            gaussian_wasserstein_distance = Inf;
        end
        error_matrix(i,j) = min(cutoff, gaussian_wasserstein_distance).^p;
    end
end
[error_location, min_idx] = min(error_matrix, [], 1);
total_local_e = sum(error_location) / X_card_larger;
opt_map = X1_labels(min_idx);

end