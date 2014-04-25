function C(version)
% correlation changes
% Non-diagonal quadratic with noise.
% Correlated

clc; close all; format long; randn('state',0); rand('state',0); 
switch version
    case 00, n = 100; v = 0.0;
    case 01, n = 100; v = 0.0;
    case 02, n = 200; v = 0.0;
    case 03, n = 300; v = 0.0;
    case 04, n = 400; v = 0.0;
    case 05, n = 500; v = 0.0;
    case 06, n = 600; v = 0.0;
    case 07, n = 700; v = 0.0;
    case 08, n = 800; v = 0.0;
    case 09, n = 900; v = 0.0;
    case 10, n =1000; v = 0.0;
        
    case 11, n = 100; v = 0.2;
    case 12, n = 200; v = 0.2;
    case 13, n = 300; v = 0.2;
    case 14, n = 400; v = 0.2;
    case 15, n = 500; v = 0.2;
    case 16, n = 600; v = 0.2;
    case 17, n = 700; v = 0.2;
    case 18, n = 800; v = 0.2;
    case 19, n = 900; v = 0.2;
    case 20, n =1000; v = 0.2;
        
    case 21, n = 100; v = 0.5;
    case 22, n = 200; v = 0.5;
    case 23, n = 300; v = 0.5;
    case 24, n = 400; v = 0.5;
    case 25, n = 500; v = 0.5;
    case 26, n = 600; v = 0.5;
    case 27, n = 700; v = 0.5;
    case 28, n = 800; v = 0.5;
    case 29, n = 900; v = 0.5;
    case 30, n =1000; v = 0.5;
  
    case 31, n = 100; v = 0.9;
    case 32, n = 200; v = 0.9;
    case 33, n = 300; v = 0.9;
    case 34, n = 400; v = 0.9;
    case 35, n = 500; v = 0.9;
    case 36, n = 600; v = 0.9;
    case 37, n = 700; v = 0.9;
    case 38, n = 800; v = 0.9;
    case 39, n = 900; v = 0.9;
    case 40, n =1000; v = 0.9;
    otherwise, return
end

p = 128; k = 5; 
lambda = 0.5*sqrt(1/n)*log(n*p); % MODIFY

alpha = 0.5;

maxit = 20; tol = 10^-6;
nrun = 100; % MODIFY
J = ones(p,nrun)==0; Ln = zeros(p,nrun); 


for run = 1:nrun
    Q = 0.5*eye(k);
    B0 = rand(k,k) < alpha;
    for i=1:k
        for j=1:i
            if B0(i,j) && i~=j
                Q(i,j) = 0.5;
            end
        end
    end
    Q = Q + Q';
    
    ord = randperm(p); J(ord(1:k),run) = 1==1; 
    X = mvnrnd(zeros(1,p),toeplitz(v.^(0:p-1)),n)';
    
    y = sum(X(J(:,run),:).*(Q*X(J(:,run),:)),1)' + randn(n,1);
    % MODIFY
    [beta,h,obj,Ln(:,run)] = acdc_QP(X,y-mean(y),lambda,maxit,tol); 
    disp(['version=' num2str(version) ' run=' num2str(run)]);
end
save(['C_' num2str(version) '.mat'],'J','Ln');
return



%% Reading the result:
clc; clear all; close all; 
prob = zeros(1,40); 
epsil = 10^-7;
for version = 1:40
    load(['C_' num2str(version) '.mat']);
    nrun = size(Ln,2); suc = 0; 
    for run = 1:nrun
        if max(Ln(~J(:,run),run)) < epsil && min(Ln(J(:,run),run)) > epsil
            suc = suc + 1;
        end
    end
    prob(version) = suc/nrun;
end
figure(2); set(gca,'FontSize',12); 
plot(100:100:1000,prob(1:10),'r.-',100:100:1000,prob(11:20),'b.-',...
    100:100:1000,prob(21:30),'g.-',100:100:1000,prob(31:40),'k.-','LineWidth',2);
legend('v=0.0','v=0.2','v=0.5','v=0.9');
xlabel('Number of Samples'); ylabel('Probability of Recovery');
title('Probability of Recovery');