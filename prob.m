function p = prob(point, mu, sigma, lambda)
%{
%% calculate probability / likelihood
% Probability that given a set of parameters \theta for the PDFs the data X can be
% observed.; this is equivalent of the likelihood of the parameters \theta
% given data points X
%
% point = [x y]
% mu = 1x2
% sigma = 2x2
% lambda = 1x2
% sigma has to be square (actually semi definite positive)
% TODO: check sigma for that and fix it if it is not.
%}
p = lambda;
for ii = 1:length(point)
    p = p * exp(-0.5 * ((point(ii,1) - mu(ii,1))./sigma(ii,ii)).^2) ./ (sqrt(2*pi) .* sigma(ii,ii));
end
end