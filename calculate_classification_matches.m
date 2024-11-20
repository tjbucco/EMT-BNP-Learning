function [match_percentage] = calculate_classification_matches(X1, X2, X1_card, X2_card, X_card_larger, n1active, n2active, C1, C2)

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
        error_matrix(i,j) = euclidean_distance;
    end
end
[error_location, min_idx] = min(error_matrix(:,:), [], 1);
opt_map = X1_labels(min_idx);
matches = 0;
if X2_card > 0
    for i = X2_card
        matches = matches + (C1(opt_map(i)) == C2(X2_labels(i)));
    end
end
match_percentage = matches./ X_card_larger;
end