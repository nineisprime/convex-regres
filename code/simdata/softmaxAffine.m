% INPUT:
%   X must be (n--by--p)
%   K is number of pieces
%   supp is the sparsity pattern

% each linear piece (beta_stars(:, k)) 
% is normalized to L2 = 2, each coordinate has some basic 
% strength

function [Y] = softmaxAffine(X, K, supp)

[n,p] = size(X);

beta_stars = randn(p,K);
beta_stars = sign(beta_stars).*(abs(beta_stars)+0.2);

%normalization
beta_stars = 2*beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));

tmp = zeros(p,K);
tmp(supp,:) = beta_stars(supp,:);
beta_stars = tmp;
%beta_stars((s+1):p,:) = 0; %only the first s coordinates are relevant

alpha_stars = randn(K,1)*0.2;

Y = softmax(X*beta_stars + ones(n,1)*alpha_stars');


% INPUT:
%    in_mat (n--by--K) matrix
%    
% OUTPUT:
%    out  (n--vec)

function [out] = softmax(in_mat)
    
out = log(sum(exp(in_mat), 2));
