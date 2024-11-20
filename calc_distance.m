function d = calc_distance(mu, sigma, lambda, mu_, sigma_, lambda_)
%{ 
This function calculates the distance between two sets of parameters. 

Input: 
    Param : old parameters
    Param_: new parameters

Output: 
    d: semi-Euclidean distance
%}
number_of_clusters = size(mu, 2);
d = 0;
for k = 1:number_of_clusters
    d = norm(mu(:, k) - mu_(:,k)) + d;
end
end