% variable selection

if (~exist('variable_selection_manual'))
    p = 50; K = 5; s = 4; n = 400; sigma=0.2;
    lambda = 0.5; mu = 0.05; maxiter = 1200;
end

beta_stars = randn(p,K);
beta_stars = sign(beta_stars).*(abs(beta_stars)+0.2);
%normalization
beta_stars = beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));

beta_stars((s+1):p,:) = 0; %only the first s coordinates are relevant

alpha_stars = randn(K,1);
alpha_stars = alpha_stars/norm(alpha_stars);

Xs = randn(n,p);
% Xs = Xs./max(max(abs(Xs)))*3;

Ys = max(Xs*beta_stars + ones(n,1)*alpha_stars', [], 2);
% Xs*-beta_stars is (n--by--K)

noise = randn(n,1)*sigma;
noisy_Ys = Ys + noise;

mean_Y = mean(noisy_Ys);
noisy_Ys = noisy_Ys - mean_Y;

lambda = lambda*sqrt(n);
[beta,z,obj] = Piecewise_Infinity(Xs',noisy_Ys, lambda, mu, maxiter);



% beta is (p--by--n)
%figure;
%plot(1:p, sum(abs(beta),2),'.r');
%[sum(beta.^2,2),sum(beta_stars.^2,2)]