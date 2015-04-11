% INPUT:
%       Sigma (p--by--p)
% 

% OUTPUT: 
%       X should be (p--by--n)
%
% mixture of uniform [-3, 3] and max-divided Gaussian
% Gaussian is truncated at 2.4
% then divided by maximum to sit between [-2.7, 2.7]

function [X, scale] = simulateBoundedGauss(p, n, unif_weight, Sigma)

X_unif = 6*rand(n,p) - 3;

if (~exist('Sigma', 'var'))
    Sigma = randn(p,p);
    Sigma = Sigma'*Sigma;

    D = diag(diag(Sigma));
    Sigma = D^(-1/2)*Sigma*D^(-1/2);
end

no_overflow = zeros(n,1);
while (sum(no_overflow) < n)
    Z = randn(10*n,p);
    X_gauss = Z*Sigma^(1/2);
    
    no_overflow = max(abs(X_gauss),[],2) < 2.7;
    disp('new attempt ...');
    num_valid = sum(no_overflow)
end

X_gauss = X_gauss(no_overflow, :);
X_gauss = X_gauss(1:n, :);

scale = 2.2/max(max(abs(X_gauss)));
X_gauss = X_gauss*scale;

mask = rand(n,1) < unif_weight;
mask = mask*ones(1,p);

X = X_unif.*mask + X_gauss.*(1-mask);



