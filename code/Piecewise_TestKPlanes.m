clc; clear all; close all; %randn('state',1); rand('state',1);
p = 100; K0 = 4; s = 4; n = 200; noise=0.1;
beta_stars = randn(p,K0); ord = randperm(p); beta_stars(ord((s+1):p),:) = 0;
alpha_stars = randn(K0,1);
Xs = randn(n,p); [Ys,c0] = max(Xs*beta_stars + ones(n,1)*alpha_stars', [], 2);
[c0,ndx] = sort(c0); Xs = Xs(ndx,:); Ys = Ys(ndx) + noise*randn(n,1);
% K-means
K = 8; c = randsample(K,n,true); lambda = max(noise*sqrt(s*log(p)/n)*n,0.01*n); mu = 1; maxit = 1000;

%Xs = 10*Xs./(ones(n,1)*sqrt(sum(Xs.^2,1)));
%Ys = sum(Xs(:,1:s),2).^2;

beta = zeros(p,K);
alpha = zeros(1,K);

for iter = 1:20
    [beta,alpha,obj,h,c,Ln] = Piecewise_KPlanesADMM(Xs',Ys,c,lambda,mu,maxit,beta,alpha);
    
    disp(['sparsity level: ' num2str(sum(max(abs(beta),[],2)<1e-5))]);
    
    disp(['Iteration ' num2str(iter) ' Objective ' num2str(obj(end))]);
    figure(1); subplot(3,2,1); plot(max(abs(beta_stars),[],2),'b*'); subplot(3,2,2); plot(Ln,'r.');
    subplot(3,2,3); plot(c0,'b.'); subplot(3,2,4); plot(c,'r.');
    subplot(3,2,[5,6]); plot(1:n,Ys,'ko',1:n,h,'r+'); drawnow;
end

% Convex
%lambda = 1; mu = 0.01;
%[beta,h,obj,Ln] = Piecewise_InfinityADMM(Xs',Ys,lambda,mu,maxit);
%figure(2); subplot(3,2,1); plot(max(abs(beta_stars),[],2),'b*'); subplot(3,2,2); plot(Ln,'r.');
%subplot(3,2,3); plot(c0,'b.'); % subplot(3,2,4); plot(c,'r.');
%subplot(3,2,[5,6]); plot(1:n,Ys,'ko',1:n,h,'r+'); drawnow;