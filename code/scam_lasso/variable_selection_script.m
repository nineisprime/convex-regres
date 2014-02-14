% variable selection

if (~exist('variable_selection_manual'))
    p = 500; K = 20; s = 10; n = 600; sigma=0.5;
    lambda = 0.3; mu = 0.05; maxiter = 1200;
end

beta_stars = randn(p,K);
beta_stars = sign(beta_stars).*(abs(beta_stars)+0.5);
%normalization
%beta_stars = beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));
if (s < p)
    beta_stars((s+1):p,:) = 0; %only the first s coordinates are relevant
end
    
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

%lambda = lambda*sqrt(n);
tic 
[beta,z,obj, Ln] = scamLasso(Xs',noisy_Ys, lambda, maxiter, 1e-5);
toc

% beta is (p--by--n)
%figure;
%plot(1:p, sum(abs(beta),2),'.r');
%[sum(beta.^2,2),sum(beta_stars.^2,2)]