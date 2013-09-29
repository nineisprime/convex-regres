clc; clear all; close all; randn('state',123); rand('state',0);
noise = 0.2; 
n = 100; p = 80; s = 3; lambda = 2*noise*sqrt(s*log(p)/n)*n; maxit = 10000; mu = 0.01;
lambda = max(lambda,1e-3*n);
prob = zeros(p,1); 
%ord = randperm(p); 
ord = 1:p;
prob(ord(1:s)) = gamrnd(10*ones(s,1),1);
prob(ord(1:s)) = prob(ord(1:s))/sum(prob(ord(1:s)));
X = randn(p,n); 
%X = X./(ones(p,1)*sqrt(sum(X.^2,1)));
%X = 0.5*rand(p,n)+1;

y = diag(X(1:s,:)'*eye(s)*X(1:s,:)) + noise*randn(n,1);

%y = exp(X'*prob) + noise*randn(n,1);
%y = log(sum(exp(X'*prob),2)) + noise*randn(n,1);

%[beta1,h1,obj1] = Piecewise_Convex(X,y,lambda,mu,maxit); % L-2-1 norm
%[beta2,h2,obj2] = Piecewise_Infinity(X,y,lambda,mu,maxit); % L-infinity-1 norm

%gamma = 1/sqrt(n);
gamma = 1e7;

[beta1,D,h1,obj1] = Quadratic_Infinity(X,y,lambda,gamma,mu,maxit); % L-2-1 norm