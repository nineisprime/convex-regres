clc; clear all; close all; randn('state',113); rand('state',0);
noise = 0.00; 
n = 2000; p = 20; s = 4; 
%lambda = max(.5*noise*sqrt(s*log(p)/n),.0001); 
lambda = 0.001*n;
%lambda = 0;
maxit = 1000; 
%mu = 0.1; 
mu = 0.01;

prob = zeros(p,1); 
%ord = randperm(p); 
ord = 1:p;
prob(ord(1:s)) = gamrnd(10*ones(s,1),1);
prob(ord(1:s)) = prob(ord(1:s))/sum(prob(ord(1:s)));

X = rand(p,n); 
%X = X./(ones(p,1)*sqrt(sum(X.^2,1)));

%X = 0.3*rand(p,n)+0.3;
%X = sign(rand(p,n)-0.5).*X;
%X = X./(sqrt(sum(X.^2,2))*ones(1,n));        

A = rand(s,s);
A = A + .5*eye(s);
A = A/norm(A,2);
%A = (1/2)*eye(s);
y = diag(X(1:s,:)'*A*X(1:s,:)) + noise*randn(n,1);

        K = 30;
        beta_stars = randn(p,K);
        
        beta_stars = beta_stars./(ones(p,1)*sqrt(sum(beta_stars.^2,1)));
        beta_stars((s+1):p,:) = 0;
        %alpha_stars = randn(K,1);
        alpha_stars = zeros(K,1);

        
        noiseless_Y = max(X'*beta_stars + ones(n,1)*alpha_stars', [], 2);
        %y = noiseless_Y + randn(n,1)*noise;       
        
%y = exp(X'*prob) + noise*randn(n,1);
%y = X'*prob + noise*randn(n,1);
%y = log(sum(exp(X'*prob),2)) + noise*randn(n,1);

%[beta1,h1,obj1] = Piecewise_Convex(X,y,lambda/sqrt(n),mu,maxit); % L-2-1 norm
lambda2 = lambda*ones(p,1);
[beta1,h1,obj1] = Piecewise_Infinity(X,y,lambda2,mu,maxit); % L-infinity-1 norm


%lambda2(1:s) = lambda/2;
%[beta2,h2,obj2] = Piecewise_Infinity(X,y,lambda2,mu,maxit); % L-infinity-1 norm

%lambda2(1:s) = lambda/10;
%[beta3,h3,obj3] = Piecewise_Infinity(X,y,lambda2,mu,maxit); % L-infinity-1 norm

figure;
subplot(3,1,1); plot(max(abs(beta1),[],2),'r.'); title('no guidance-weight');
%subplot(3,1,2); plot(max(abs(beta3),[],2),'r.'); title('2x guidance-weight');
%subplot(3,1,3); plot(max(abs(beta3),[],2),'r.'); title('10x guidance-weight');

set(gca,'FontSize',12);

