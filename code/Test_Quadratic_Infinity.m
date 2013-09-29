clc; clear all; close all; 

%randn('state',123); rand('state',0);
noise = 0; 
n = 70; p = 20; s = 3; lambda = 2*noise*sqrt(s*log(p)/n)*n; 

maxit = 5000; mu = .1;
lambda = max(lambda,1e-3*n);

X = randn(p,n); 


%X = X./(ones(p,1)*sqrt(sum(X.^2,1)));
%X = 0.5*rand(p,n)+1;

K = 10;
beta_star = randn(p,K);
beta_star((s+1):p,:)=0;

tmp = randn(s,s);
A = tmp*tmp';
A = A/norm(A,2);
%A = eye(s);

y = diag(X(1:s,:)'*A*X(1:s,:)) + max(X'*beta_star,[],2) + noise*randn(n,1);

%y = exp(X'*prob) + noise*randn(n,1);
%y = log(sum(exp(X'*prob),2)) + noise*randn(n,1);

%[beta1,h1,obj1] = Piecewise_Convex(X,y,lambda,mu,maxit); % L-2-1 norm
%[beta2,h2,obj2] = Piecewise_Infinity(X,y,lambda,mu,maxit); % L-infinity-1 norm

%gamma = 1/sqrt(n);
%gamma = 1/sqrt(n);
%gamma = 1/n;
gamma = 1;

[beta1,D,h1,obj1] = Quadratic_Infinity(X,y,lambda,gamma,mu,maxit); % L-2-1 norm