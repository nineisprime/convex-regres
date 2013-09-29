% variable selection

if (~exist('variable_selection_manual'))
    p = 100; K = 10; s = 4; n = 700; sigma=0;
    lambda = 0.001; mu = 0.01; maxiter = 1000;
end

beta_stars = randn(p,K);
beta_stars = sign(beta_stars).*(abs(beta_stars)+0.2);
%normalization
beta_stars = beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));

beta_stars((s+1):p,:) = 0; %only the first s coordinates are relevant

alpha_stars = randn(K,1);
alpha_stars = alpha_stars/norm(alpha_stars);

Xs = randn(n,p);
Xs = Xs./max(max(abs(Xs)))*3;

Ys = max(Xs*beta_stars + ones(n,1)*alpha_stars', [], 2);
% Xs*-beta_stars is (n--by--K)

Ys = Ys + randn(n,1)*sigma;
Ys = Ys - mean(Ys);

[beta,z,obj] = SCAM_QP(Xs',Ys, lambda, maxiter, 1e-4);

% beta is (p--by--n)
%figure;
%plot(1:p, sum(abs(beta),2),'.r');
%[sum(beta.^2,2),sum(beta_stars.^2,2)]