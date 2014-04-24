function C1(version)
% Non-diagonal quadratic with noise.
clc; close all; format long; randn('state',0); rand('state',0); 
switch version
    case 00, n = 500; v = 0.8;
    case 01, n = 100; v = 0.2;
    case 02, n = 200; v = 0.2;
    case 03, n = 300; v = 0.2;
    case 04, n = 400; v = 0.2;
    case 05, n = 500; v = 0.2;
    case 06, n = 600; v = 0.2;
    case 07, n = 700; v = 0.2;
    case 08, n = 800; v = 0.2;
    case 09, n = 900; v = 0.2;
    case 10, n =1000; v = 0.2;
        
    case 11, n = 100; v = 0.4;
    case 12, n = 200; v = 0.4;
    case 13, n = 300; v = 0.4;
    case 14, n = 400; v = 0.4;
    case 15, n = 500; v = 0.4;
    case 16, n = 600; v = 0.4;
    case 17, n = 700; v = 0.4;
    case 18, n = 800; v = 0.4;
    case 19, n = 900; v = 0.4;
    case 20, n =1000; v = 0.4;
        
    case 21, n = 100; v = 0.6;
    case 22, n = 200; v = 0.6;
    case 23, n = 300; v = 0.6;
    case 24, n = 400; v = 0.6;
    case 25, n = 500; v = 0.6;
    case 26, n = 600; v = 0.6;
    case 27, n = 700; v = 0.6;
    case 28, n = 800; v = 0.6;
    case 29, n = 900; v = 0.6;
    case 30, n =1000; v = 0.6;
        
    case 31, n = 100; v = 0.8;
    case 32, n = 200; v = 0.8;
    case 33, n = 300; v = 0.8;
    case 34, n = 400; v = 0.8;
    case 35, n = 500; v = 0.8;
    case 36, n = 600; v = 0.8;
    case 37, n = 700; v = 0.8;
    case 38, n = 800; v = 0.8;
    case 39, n = 900; v = 0.8;
    case 40, n =1000; v = 0.8;
    otherwise, return
end

p = 128; k = 5; lambda = 4*sqrt(log(n*p)/n); maxit = 20; tol = 10^-6;
nrun = 200; J = ones(p,nrun)==0; Ln = zeros(p,nrun);
alpha = 0.5;
if alpha<10^-16, Q = eye(k); else Q = 0.5*eye(k); B0 = rand(k,k) < alpha;
for i = 1:k, for j = 1:i, if B0(i,j) && i~=j, Q(i,j) = 0.5; end; end; end; Q = Q + Q';end

for run = 1:nrun
    ord = randperm(p); J(ord(1:k),run) = 1==1; 
    X = mvnrnd(zeros(1,p),toeplitz(v.^(0:p-1)),n)';
    y = sum(X(J(:,run),:).*(Q*X(J(:,run),:)),1)' + randn(n,1);
    [beta,h,obj,Ln(:,run)] = SCAM_QP(X,y-mean(y),lambda,maxit,tol);
    disp(['version=' num2str(version) ' run=' num2str(run)]);
end
save(['SCAM/C1_' num2str(version) '.mat'],'J','Ln');
return