clc; clear all; close all; randn('state',123); rand('state',0);

noise = 0.005;
s = 4; K = 20;

p_ls = [50,100];
n_ls = [100,200,300,400,500,600,700,800];
succ_prob_ls = zeros(length(n_ls),length(p_ls));
maxit = 1500; mu = 0.01;

num_trials = 10;

for ii=1:length(n_ls)
for ij=1:length(p_ls)
    n = n_ls(ii);
    p = p_ls(ij);
    lambda = max(0.2*noise*sqrt(s*log(p)/n)*n, 0); 
    
    extra_weight = 0;
    for it = 1:num_trials
        disp(sprintf('ii:%d   ii:%d     extra_weight:%d\n', ii, it, extra_weight));
        succ_prob_ls'
        
        beta_stars = randn(p,K);
        beta_stars = beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));
        beta_stars((s+1):p,:) = 0;
        %alpha_stars = randn(K,1);
        alpha_stars = zeros(K,1);
        
        X = randn(n,p);
        %X = X./(ones(n,1)*sqrt(sum(X.^2,1)));
        X = X./5;
        
        Y = diag(X(:,1:s)*eye(s)*X(:,1:s)') + noise*randn(n,1);
        %noiseless_Y = max(X*beta_stars + ones(n,1)*alpha_stars', [], 2);
        %Y = noiseless_Y + randn(n,1)*noise;
        
        [beta2,h2,obj2] = Piecewise_Infinity(X',Y,lambda,mu,maxit);
        %extra_weight = extra_weight+sum(max(abs(beta2((s+1):p,n)),[],2));
        
        %prods = X*beta_stars;
        %[vals,ixs] = max(prods,[],2);
        %beta_golds = beta_stars(:,ixs);
        %extra_weight = extra_weight+sum(max(abs(beta_golds - beta2),[],2))/sum(max(abs(beta_golds),[],2));
        extra_weight = extra_weight+sum(max(abs(beta2(1:s,:)),[],2))/sum(max(abs(beta2),[],2));
        
        close all;
        %[beta1,h1,obj1] = Piecewise_Convex(X,y,lambda/sqrt(n),mu,maxit);        
    end
    succ_prob_ls(ii,ij) = extra_weight/num_trials;
    save 'tmp_recovery3.mat'
end
end