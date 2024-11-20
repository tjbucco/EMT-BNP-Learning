function [total_local_e, opt_map] = calculate_local_error(X1, X2, X1_card, X2_card, X_card_larger, n1active, n2active, cutoff, p)

error_matrix = NaN(X1_card, X2_card);
X1_labels = find(n1active);
X2_labels = find(n2active);
for i = 1:X1_card
    for j = 1:X2_card
        jointly_exist = ~isnan(X1(1,X1_labels(i))) & ~isnan(X2(1,X2_labels(j)));
        if any(jointly_exist)
            euclidean_distance = norm(X1(1:2, X1_labels(i)) - X2(1:2, X2_labels(j)));
        else
            euclidean_distance = Inf;
        end
        error_matrix(i,j) = min(cutoff, euclidean_distance).^p;
    end
end
[error_location, min_idx] = min(error_matrix, [], 1);
total_local_e = sum(error_location) / X_card_larger;
opt_map = X1_labels(min_idx);

end