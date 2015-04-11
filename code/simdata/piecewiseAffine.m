% INPUT:
%   X must be (n--by--p)
%   K is number of pieces
%   s is the sparsity level

% each linear piece (beta_stars(:, k)) 
% is normalized to L2 = 2, each coordinate has some basic 
% strength

function [Y] = piecewiseAffine(X, K, s)

[n,p] = size(X);

beta_stars = randn(p,K);
beta_stars = sign(beta_stars).*(abs(beta_stars)+0.1);

%normalization
beta_stars = 2*beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));

beta_stars((s+1):p,:) = 0; %only the first s coordinates are relevant

alpha_stars = randn(K,1)*0.3;

Y = max(X*beta_stars + ones(n,1)*alpha_stars', [], 2);

