function C(version)
% correlation changes
% Non-diagonal quadratic with noise.
% Correlated
%
% bounded Gaussian and uniform mixture
%
addpath('../simdata');

clc; close all; format long; randn('state',0); rand('state',0); 
switch version
    case 01, n = 200; v = 0.0;
    case 02, n = 400; v = 0.0;
    case 03, n = 600; v = 0.0;
    case 04, n = 800; v = 0.0;
    case 05, n = 1000; v = 0.0;
    case 06, n = 1200; v = 0.0;
    case 07, n = 1400; v = 0.0;
        
    case 08, n = 200; v = 0.2;
    case 09, n = 400; v = 0.2;
    case 10, n = 600; v = 0.2;
    case 11, n = 800; v = 0.2;
    case 12, n = 1000; v = 0.2;
    case 13, n = 1200; v = 0.2;
    case 14, n = 1400; v = 0.2;
        
    case 15, n = 200; v = 0.5;
    case 16, n = 400; v = 0.5;
    case 17, n = 600; v = 0.5;
    case 18, n = 800; v = 0.5;
    case 19, n = 1000; v = 0.5;
    case 20, n = 1200; v = 0.5;
    case 21, n = 1400; v = 0.5;
  
    case 22, n = 200; v = 0.9;
    case 23, n = 400; v = 0.9;
    case 24, n = 600; v = 0.9;
    case 25, n = 800; v = 0.9;
    case 26, n = 1000; v = 0.9;
    case 27, n = 1200; v = 0.9;
    case 28, n = 1400; v = 0.9;
    otherwise, return
end
num_versions = 28;

p = 128; k = 5; 
SNR = 5;

lambda = 0.5*sqrt(1/n)*log(n*p); % MODIFY

alpha = 0.5;

maxit = 20; tol = 10^-6;
nrun = 40; % MODIFY
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
    %X = mvnrnd(zeros(1,p),toeplitz(v.^(0:p-1)),n)';
    Sigma = toeplitz(v.^(0:p-1));
    unif_weight = 0.05;
    X = simulateBoundedGaussCopula(p, n, unif_weight, Sigma);
    X = X';
    y = sum(X(J(:,run),:).*(Q*X(J(:,run),:)),1)'; 
    y = y - mean(y);
    y = SNR*y / std(y);
    
    y = y+ randn(n,1);
    % MODIFY
    [beta,h,obj,Lnvex,Lncave] = acdc_QP(X,y-mean(y),lambda,maxit,tol); 
    
    epsil = 1e-6;
    succ1 = min(Ln(J(:,run),run)) > epsil
    succ2 = sum(Ln(:,run) > epsil) < 11    
    
    Ln(:,run) = max(Lnvex, Lncave);
    disp(['version=' num2str(version) ' run=' num2str(run)]);
end
%save(['C_' num2str(version) '.mat'],'J','Ln');
return



%% Reading the result:
clc; clear all; close all; 
prob = zeros(1,num_versions); 
epsil = 1e-6;
for version = 1:num_versions
    load(['C_' num2str(version) '.mat']);
    nrun = size(Ln,2); suc = 0; 
    for run = 1:nrun
        if min(Ln(J(:,run), run)) > epsil && ...
           sum(Ln(:,run) > epsil) < 11
            suc = suc + 1;
        end
    end
    prob(version) = suc/nrun;
end

% much of the figures still need to be edited

figure(2); set(gca,'FontSize',12); 
plot(100:100:1000,prob(1:10),'r.-',100:100:1000,prob(11:20),'b.-',...
    100:100:1000,prob(21:30),'g.-',100:100:1000,prob(31:40),'k.-','LineWidth',2);

legend('v=0.0','v=0.2','v=0.5','v=0.9');
xlabel('Number of Samples'); 
ylabel('Probability of Recovery');
title('Probability of Recovery');
%set(gca, 'LooseInset', get(gca, 'TightInset'));

tightInset = get(gca, 'TightInset');
position(1) = tightInset(1);
position(2) = tightInset(2);
position(3) = 1 - tightInset(1) - tightInset(3);
position(4) = 1 - tightInset(2) - tightInset(4);
set(gca, 'Position', position);
set(gca,'units','centimeters')
pos = get(gca,'Position');

ti = get(gca,'TightInset');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

h = gcf;
%saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/Curve4.pdf');