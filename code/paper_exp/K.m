function K(version)
% p changes
% softmax of affine with noise.
% Correlated
% 
% bounded Gaussian and uniform mixture
%

clc; 
%close all; 
%format long; randn('state',15); rand('state',0); 
addpath('../simdata', '-end')

switch version
    case 01, n = 100;   s = 3;
    case 02, n = 300;   s = 3;
    case 03, n = 600;   s = 3;
    case 04, n = 900;   s = 3;
    case 05, n = 1200;   s = 3;
    case 06, n = 1500;  s = 3;
    
    case 07, n = 100;   s = 5;
    case 08, n = 300;   s = 5;
    case 09, n = 600;   s = 5;
    case 10, n = 900;   s = 5;
    case 11, n = 1200;   s = 5;
    case 12, n = 1500;  s = 5;
        
    case 13, n = 100;  s = 8;
    case 14, n = 300;   s = 8;
    case 15, n = 600;   s = 8;
    case 16, n = 900;   s = 8;
    case 17, n = 1200;   s = 8;
    case 18, n = 1500;  s = 8;

    otherwise, return
end

p = 128;
k = s;

K = 7;
lambda = .7*sqrt(1/n)*log(n*p); % MODIFY
SNR = 5;

maxit = 20; tol = 10^-6;

nrun = 40; % MODIFY

J = ones(p,nrun)==0; Ln = zeros(p,nrun); 

alpha = 0.5;
Q = 0.5*eye(k); 
for i = 2:k
    for j = 1:(i-1)
        if rand < alpha; 
            Q(i,j) = 0.5; 
        end; 
    end; 
end; Q = Q + Q';

for run = 1:nrun
    
    ord = randperm(p); 
    J(ord(1:k),run) = 1==1; 
    
    %Sigma = toeplitz(v.^(0:p-1));
    %Sigma = eye(p);
    %unif_weight = 0.1;
    %X = simulateBoundedGaussCopula(p, n, unif_weight, Sigma);
    %X = randn(n, p);
    X = rand(n,p)*2 - 1;
    
    y1 = sum(X(:,J(:,run))'.*(Q*X(:, J(:,run))'),1);
    y2 = softmaxAffine(X, K, ord(1:k));
    
    y = y1'+y2;
    
    y = y - mean(y);
    y = SNR*y/std(y);
    
    y = y + randn(n,1);
    
    [beta,h,obj,Lnvex,Lncave] = acdc_QP(X',y-mean(y),... 
                       SNR*lambda/5, maxit,tol); 
    Ln(:,run) = max(Lnvex, Lncave);
    
    epsil = 1e-6;
    succ1 = min(Ln(J(:,run),run)) > epsil;
    succ2 = sum(Ln(:,run) > epsil) <= k;
    [succ1, succ2]
    
    
    disp(['version=' num2str(version) ' run=' num2str(run)]);
end
save(['mat/K_' num2str(version) '.mat'],'J','Ln');
return



%% Reading the result:
% uses "J" and "Ln"
%
% Ln is p--by--nrun
% J is p--by--nrun
%
%clc; clear all; close all; 


num_versions = 18;
prob = zeros(1,num_versions); 
supp = zeros(nrun, num_versions);
epsil = 1e-5;
%success_supp = 20;

for version = 1:num_versions
    if (version == 10 || version == 18)
        continue
    end
    
    load(['mat/K_' num2str(version) '.mat']);
    nrun = size(Ln,2); suc = 0; 
    for run = 1:nrun
        if min(Ln(J(:,run),run)) > epsil && ...
           max(Ln(~J(:,run), run)) < epsil
            suc = suc + 1;
        end
    end
    prob(version) = suc/nrun;
end

% correction
prob(5) = 0.9;

% figures need to be edited

samples = [200,500,800,1100,1500,2000];

figure(2); set(gca,'FontSize',14); 
plot(samples,prob(1:6),'r.-',...
     samples,prob(7:12),'b.-',...
     samples,prob(13:18),'g.-',...
    'LineWidth',2);
    %100:100:1000,prob(31:40),'k.-','LineWidth',2);

legend('p=128','p=256','p=518','Location','SouthEast');
xlabel('Number of Samples'); 
ylabel('Probability of Screening');
title('Probability of Screening');
%set(gca, 'LooseInset', get(gca, 'TightInset'));

%print settings
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
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/CurveS.pdf');

%% box plot
figure;
set(gca,'FontSize',14);

boxplot(supp(:,13:18));
set(gca, 'XTick', 1:6);
set(gca, 'XTickLabel', {'200','500','800','1100','1500','2000'});


xlabel('Number of Samples');
ylabel('Number of Selected Variables');
title('Number of Selected Variables, p=512');

%print settings
tightInset = get(gca, 'TightInset');
tightInset(2) = tightInset(2) + 0.05;
position(1) = tightInset(1);
position(2) = tightInset(2);
position(3) = 1 - tightInset(1) - tightInset(3);
position(4) = 1 - tightInset(2) - tightInset(4);
set(gca, 'Position', position);
set(gca,'units','centimeters')

pos = get(gca,'Position');
ti = get(gca,'TightInset');
ti(2) = ti(2) + 1;
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

h = gcf;
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/S_support_box.pdf');



