% INPUT:
%       Sigma (p--by--p)
% 

% OUTPUT: 
%       X should be (p--by--n)
%
% mixture of uniform [-3, 3] and a mixture of Gaussians
%
% means are drawn from N(0, Sigma)*1/2
% std of each component is uniform(0.3, 0.8)
% weights are drawn from uniform distribution

function [X] = simulateBoundedGaussMixture(p, n, unif_weight, Sigma, T)

border = 4;
X_unif = 2*border*rand(n,p) - border;
mix_border = 3.7;

% generate the component means
while (1)
    means = 0.5*randn(T,p);
    means = means*Sigma^(1/2);
    
    
    if (max(max(abs(means))) < 3)
        means = 2.5*means/max(max(abs(means)));
        means
        break;
    end
end

% generate the component weights
weights = rand(T);
weights = weights/sum(weights);
weights = floor(weights*n);
weights(1) = n - sum(weights(2:T));

weights

X_mix = zeros(0, p);

for (ii = 1:T)
    
    cur_sd = rand*0.5 + 0.3;
    
    while (1)
        X_comp = randn(2*n, p)*cur_sd + ones(2*n, 1)*means(ii,:);
        no_overflow = max(abs(X_comp), [], 2) < mix_border;
        if (sum(no_overflow) > weights(ii))
            break;
        end
    end
    
    X_comp = X_comp(no_overflow,:);
    X_comp = X_comp(1:weights(ii), :);
    X_mix = [X_comp; X_mix];
end

mask = rand(n,1) < unif_weight;
mask = mask*ones(1,p);

X = X_unif.*mask + X_mix.*(1-mask);



