%=====================================%
% Sparse Convex Additive Model
% Author name withheld
% Last Edit: 5/31/2013
%
% use the MOSEK matlab toolbox (academic license availabe for free) for best effect
% uses the MOSEK quadprog function

% main function is SCAM_QP.m
% T.m file contains code for synthetic experiment

% Example usage:

n = 600;
p = 100;
k = 5; lambda = 4*sqrt(log(n*p)/n); maxit = 20; tol = 10^-6; Q = eye(k);
X = randn(p,n); y = sum(X(1:k,:).*(Q*X(1:k,:)),1)' + randn(n,1);
[beta,h,obj] = SCAM_QP(X,y-mean(y),lambda,maxit,tol);

% magnitudes is p-dimensional, large magnitudes == relevant variables
magnitudes = max(abs(beta),[],2);
