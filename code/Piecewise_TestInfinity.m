%clc; clear all; close all; 
%randn('state',123); rand('state',0);

addpath('simdata');

noise = 0.2; 
n = 1000; 
p = 30; 
s = 3; 

lambda = noise*sqrt(s*log(p)/n)*n; 

maxit = 2000; 
mu = 0.005;

lambda = max(lambda,1e-3*n);

X = 2*rand(p,n) - 1;

K = 4;
y = softmaxAffine(X', K, 1:s) + noise*randn(n,1);


%[beta1,h1,obj1] = Piecewise_Convex(X,y,lambda,mu,maxit); % L-2-1 norm
[beta2,h2,obj2] = Piecewise_Infinity(X,y,lambda,mu,maxit); % L-infinity-1 norm



%gamma = 1/sqrt(n);
%gamma = 1e7;


%[beta1,D,h1,obj1] = Quadratic_Infinity(X,y,lambda,gamma,mu,maxit); % L-2-1 norm


