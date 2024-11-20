function [mu, sigma, lambda] = maximization(Data, mu, sigma, lambda)
%{ 
This function calculates the second step of the EM algorithm, Maximization.
It updates the parameters of the Normal distributions according to the new 
labled dataset.

Input: 
    Data : nx3 (number of data points , [x, y, label])
    Param: (mu, sigma, lambda)

Output: 
    Param: updated parameters 
%}

points_in_cluster = histc(Data(end, :), 1:size(mu, 2));

% calculate the weights
lambda = points_in_cluster./sum(points_in_cluster);

for k = find(lambda > 0)

    % calculate the means
    mu(:, k) = mean(Data(1:end - 1, Data(end, :) == k), 2);

    % calculate the variances
%     diag_elem = zeros(size(Data, 1) - 1, 1);
%     for j = 1:size(Data, 1) - 1
%         diag_elem(j) = std(Data(j, Data(end,:) == k));
%     end
%     if sum(diag_elem) > 0
%         sigma(:,:,k) = diag(diag_elem);
%     end
end

end