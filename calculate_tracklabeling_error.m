function error = calculate_tracklabeling_error(X1, X2, X1_card, X2_card, n1active, n2active, opt_map, assignment)

error = 0;
X1_labels = find(n1active);
X2_labels = find(n2active);

for i = 1:X2_card
    error = 1 - assignment(X2_labels(i), opt_map(i)) + error;
end

end