% INPUT:
%       Sigma (p--by--p)
% 

% OUTPUT: 
%       X should be (p--by--n)
%
% mixture of uniform [-2, 2] and max-divided Gaussian
% Gaussian is truncated at 2.4
% then divided by maximum to sit between [-2.7, 2.7]

function [X] = simulateBoundedGaussCopula(p, n, unif_weight, Sigma)

border = 2;
mix_border = 1.8;

X_unif = 2*border*rand(n,p) - border;

if (~exist('Sigma', 'var'))
    Sigma = randn(p,p);
    Sigma = Sigma'*Sigma;

    D = diag(diag(Sigma));
    Sigma = D^(-1/2)*Sigma*D^(-1/2);
end


Z = randn(n,p);
X_gauss = Z*Sigma^(1/2);
    
X_gauss = normcdf(X_gauss);
a = 5;
X_gauss = 2*mix_border*(exp(a/2) + 1)/(exp(a/2)-1) * ...
          (exp(a*(X_gauss-1/2))./(exp(a*(X_gauss-1/2)) + 1) ...
          - 1/(exp(a/2)+1)) ...
          - mix_border;


mask = rand(n,1) < unif_weight;
mask = mask*ones(1,p);

X = X_unif.*mask + X_gauss.*(1-mask);


%a = 5;

%x = 0:0.01:1;
%fx = (exp(a/2) + 1)/(exp(a/2)-1) * ...
%          (exp(a*(x-1/2))./(exp(a*(x-1/2)) + 1) ...
%          - 1/(exp(a/2)+1));
%plot(fx)