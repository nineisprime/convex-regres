function S(version)
% p changes
% softmax of affine with noise.
% Correlated
% 
% bounded Gaussian and uniform mixture
%

clc; 
%close all; 
format long; randn('state',15); rand('state',0); 
addpath('../simdata', '-end')

switch version
    case 01, n = 400; p = 128;
    case 02, n = 800; p = 128;
    case 03, n = 1200; p = 128;
    case 04, n = 1600; p = 128;
    case 05, n = 2000; p = 128;
    case 06, n = 2400; p = 128;
        
    case 07, n = 400; p = 256;
    case 08, n = 800; p = 256;
    case 09, n = 1200; p = 256;
    case 10, n = 1600; p = 256;
    case 11, n = 2000; p = 256;
    case 12, n = 2400; p = 256;
        
    case 13, n = 400; p = 512;
    case 14, n = 800; p = 512;
    case 15, n = 1200; p = 512;
    case 16, n = 1600; p = 512;
    case 17, n = 2000; p = 512;
    case 18, n = 2400; p = 512;
        
    otherwise, return
end

v = 0.5;
k = 5; 
K = 7;
lambda = 0.5*sqrt(1/n)*log(n*p); % MODIFY
SNR = 5;

maxit = 20; tol = 10^-6;

nrun = 20; % MODIFY

J = ones(p,nrun)==0; Ln = zeros(p,nrun); 


for run = 1:nrun
    
    ord = randperm(p); 
    J(ord(1:k),run) = 1==1; 
    
    Sigma = toeplitz(v.^(0:p-1));
    unif_weight = 0.1;
    X = simulateBoundedGaussCopula(p, n, unif_weight, Sigma);
    
    %y = sum(X(J(:,run),:).*(Q*X(J(:,run),:)),1)' + 
    %randn(n,1);
    y = softmaxAffine(X, K, ord(1:k));
    y = y - mean(y);
    y = SNR*y/std(y);
    
    y = y + randn(n,1);
    
    [beta,h,obj,Lnvex,Lncave] = acdc_QP(X',y-mean(y),... 
                       SNR*lambda/5, maxit,tol); 
    Ln(:,run) = max(Lnvex, Lncave);
    
    epsil = 1e-6;
    succ1 = min(Ln(J(:,run),run)) > epsil
    succ2 = sum(Ln(:,run) > epsil) < 11
    
    
    disp(['version=' num2str(version) ' run=' num2str(run)]);
end
save(['mat/S_' num2str(version) '.mat'],'J','Ln');
return



%% Reading the result:
%clc; clear all; close all; 
num_versions = 18;
prob = zeros(1,num_versions); 
epsil = 1e-6;

for version = 1:num_versions
    load(['mat/S_' num2str(version) '.mat']);
    nrun = size(Ln,2); suc = 0; 
    for run = 1:nrun
        if sum(Ln(:,run) > epsil) < 11 && ... 
           min(Ln(J(:,run),run)) > epsil
            suc = suc + 1;
        end
    end
    prob(version) = suc/nrun;
end

figure(2); set(gca,'FontSize',12); 
plot(100:100:1000,prob(1:5),'r.-',...
    100:100:1000,prob(5:10),'b.-',...
    100:100:1000,prob(11:15),'g.-',...
    'LineWidth',2);
    %100:100:1000,prob(31:40),'k.-','LineWidth',2);

legend('p=128','p=256','p=518');
xlabel('Number of Samples'); 
ylabel('Probability of Screening');
title('Probability of Screening');
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


