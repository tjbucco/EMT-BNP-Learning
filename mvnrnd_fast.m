function [r] = mvnrnd_fast(mu,sigma)
%MVNRND Random vectors from the multivariate normal distribution.
%   R = MVNRND(MU,SIGMA) returns an N-by-D matrix R of random vectors
%   chosen from the multivariate normal distribution with mean vector MU,
%   and covariance matrix SIGMA.  MU is an N-by-D matrix, and MVNRND
%   generates each row of R using the corresponding row of MU.  SIGMA is a
%   D-by-D symmetric positive semi-definite matrix, or a D-by-D-by-N array.
%   If SIGMA is an array, MVNRND generates each row of R using the
%   corresponding page of SIGMA, i.e., MVNRND computes R(I,:) using MU(I,:)
%   and SIGMA(:,:,I).  If the covariance matrix is diagonal, containing
%   variances along the diagonal and zero covariances off the diagonal,
%   SIGMA may also be specified as a 1-by-D matrix or a 1-by-D-by-N array,
%   containing just the diagonal. If MU is a 1-by-D vector, MVNRND
%   replicates it to match the trailing dimension of SIGMA.
%
%   R = MVNRND(MU,SIGMA,N) returns a N-by-D matrix R of random vectors
%   chosen from the multivariate normal distribution with 1-by-D mean
%   vector MU, and D-by-D covariance matrix SIGMA.
%
%   Example:
%
%      mu = [1 -1]; Sigma = [.9 .4; .4 .3];
%      r = mvnrnd(mu, Sigma, 500);
%      plot(r(:,1),r(:,2),'.');
%
%   See also MVTRND, MVNPDF, MVNCDF, NORMRND.

%   R = MVNRND(MU,SIGMA,N,T) supplies the Cholesky factor T of
%   SIGMA, so that SIGMA(:,:,J) == T(:,:,J)'*T(:,:,J) if SIGMA is a 3D array or SIGMA
%   == T'*T if SIGMA is a matrix.  No error checking is done on T.
%
%   [R,T] = MVNRND(...) returns the Cholesky factor T, so it can be
%   re-used to make later calls more efficient, although there are greater
%   efficiency gains when SIGMA can be specified as a diagonal instead.

%   Copyright 1993-2020 The MathWorks, Inc.

[n,d] = size(mu);

% Just the diagonal of sigma has been specified.
t = sqrt(sigma);
r = randn(n,d).*t + mu;