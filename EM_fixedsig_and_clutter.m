function [Data_f, mu_f, sigma_f, lambda_f] = EM_fixedsig_and_clutter(Data, mu, sigma, lambda)
%{ 
This is the main EM algorithm. It has two steps Expectation and
Maximization. The whole process is done in a while loop until a desired
error has reached. 

Input: 
    Data: the dataset including the points and labels [x y label]
    mu, sigma, lambda: parameters of the Normal Distributions

Output: 
    Data_f: the dataset with updated labels
    mu_f, sigma_f, lambda_f: final parameters that the algorithm has converged to
%}

shift = 10000;  % a big number
iter = 0;       % counter
epsilon = 0.001; % percision
formatSpec = 'iteration: %d, error: %2.4f \n';
%mu_f, sigma_f, lambda_f

while shift > epsilon
    iter = iter + 1;

    % E-step
    Data_ = expectation(Data, mu, sigma, lambda);
    
    % M-step
    [mu_, ~, lambda_] = maximization(Data_, mu, sigma, lambda);
    mu_(:, end) = mu(:,end);
    
    % calculate the distance/error from the previous set of params
    shift = calc_distance(mu, sigma, lambda, mu_, sigma, lambda_);
    
    %fprintf(formatSpec, iter, shift);    

    Data = Data_;
    mu = mu_;
    lambda = lambda_;

    clear Data_ mu_ lambda_
end
Data_f = Data;
mu_f = mu;
sigma_f = sigma;
lambda_f = lambda;
end