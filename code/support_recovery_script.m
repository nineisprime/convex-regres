%clc; clear all; close all; randn('state',123); rand('state',0);

noise = 0.05;
s = 3; K = 6;

lambda0 = 1;
%lambda_ls = [lambda0/30, lambda0/10, lambda0/3, lambda0, lambda0*3];
lambda_ls = 0.5;

p_ls = 50;
n_ls = 1200:100:1400;
succ_prob_ls = zeros(length(n_ls),length(p_ls),length(lambda_ls));
false_discovery = zeros(length(n_ls),length(p_ls),length(lambda_ls));
false_neglect = zeros(length(n_ls),length(p_ls),length(lambda_ls));

maxit = 2000; mu = 0.1;

num_trials = 5;

for ii=1:length(n_ls)
for ij=1:length(p_ls)
    n = n_ls(ii);
    p = p_ls(ij);
   
    %lambda0 = max(noise*sqrt(s*log(p)/n)*n, 1e-5); 
   
    
    
    for il=1:length(lambda_ls)
    num_succ=0;
    num_false_discovery=0;
    num_false_neglect=0;
    for it = 1:num_trials
        disp(sprintf('ii:%d   it:%d     num_succ:%d\n', ii, it, num_succ));
        
        lambda = lambda_ls(il);
        beta_stars = randn(p,K);
        beta_stars = beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));
        beta_stars((s+1):p,:) = 0;
        %alpha_stars = randn(K,1);
        alpha_stars = zeros(K,1);
        
        X = 2*rand(n,p)-1;
        %X = X./(sqrt(sum(X.^2,2))*ones(1,p));
        
        y_lin = max(X*beta_stars + ones(n,1)*alpha_stars', [], 2);
        
        A = randn(s,s);
        A = A*A' + .5*eye(s);
        A = A/norm(A,2);
        y_quad = diag(X(:,1:s)*A*X(:,1:s)');
        
        y_cube = sum(abs(X(:,1:(s-1))).^3,2);
        
        Y = .5*y_cube+.3*y_lin + .5*y_quad + randn(n,1)*noise;
        
        %lambda = lambda/10;
        [beta2,h2,obj2] = Piecewise_Infinity(X',Y,lambda,mu,maxit);
        lambda2 = ones(p,1)*lambda ./ ((max(abs(beta2),[],2) + 1e-6*ones(p,1)));
        %lambda2 = lambda2/min(lambda2)*lambda/2;
        
        for jj=1:4
            if (jj > 3)
                c=2.5;
            else
                c=1;
            end
            [beta2,h2,obj2] = Piecewise_Infinity(X',Y,lambda2,mu,c*maxit);
            lambda2 = ones(p,1)*lambda ./ ((max(abs(beta2),[],2) + 1e-6*ones(p,1))); 
            %lambda2 = lambda3/min(lambda3)*lambda/2;
        end
        
        no_false_positive = max(max(abs(beta2((s+1):p,:)),[],2)) < 1e-4;
        no_false_negative = min(max(abs(beta2(1:s,:)),[],2)) > 1e-4;
        if (no_false_positive && no_false_negative)
            num_succ = num_succ+1;
            num_false_discovery = num_false_discovery+1;
            num_false_neglect = num_false_neglect+1;
        end
        close all;
        %[beta1,h1,obj1] = Piecewise_Convex(X,y,lambda/sqrt(n),mu,maxit);        
    end
    succ_prob_ls(ii,ij,il) = num_succ/num_trials;
    false_discovery(ii,ij,il) = num_false_discovery;
    false_neglect(ii,ij,il) = num_false_neglect;
    succ_prob_ls
    %save 'tmp_recovery_reweight3.mat'
    end
    

   
end
end