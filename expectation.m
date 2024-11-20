function Data = expectation(Data, mu, sigma, lambda)
%{
This function calculates the first step of the EM algorithm, Expectation.
It calculates the probability of each specific data point belong to each
cluster or class

Input: 
    Data : nx3 (number of data points , [x, y, label])
    Param: (mu, sigma, lambda)

Output: 
    Data: the dataset with updated labels
%}
number_of_clusters = size(mu, 2);
p_cluster = NaN(1, number_of_clusters);
for ii = 1: size(Data,2)
    x = Data(1:end - 1, ii);
    for k = 1:size(mu, 2)
        p_cluster(1, k) = prob(x, mu(:, k), sigma(:,:,k), lambda(1,k));
    end
    [~, Data(end, ii)] = max(p_cluster);
    %     if p_cluster1 > p_cluster2
    %         Data(ii, 3) = 1;
    %     else
    %         Data(ii, 3) = 2;
    %     end
end
end