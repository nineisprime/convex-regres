function S(version)
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
    case 01, n = 200; p = 128;
    case 02, n = 500; p = 128;
    case 03, n = 800; p = 128;
    case 04, n = 1100; p = 128;
    case 05, n = 1500; p = 128;
    case 06, n = 2000; p = 128;
        
    case 07, n = 200; p = 256;
    case 08, n = 500; p = 256;
    case 09, n = 800; p = 256;
    case 10, n = 1100; p = 256;
    case 11, n = 1500; p = 256;
    case 12, n = 2000; p = 256;
        
    case 13, n = 200; p = 512;
    case 14, n = 500; p = 512;
    case 15, n = 800; p = 512;
    case 16, n = 1100; p = 512;
    case 17, n = 1500; p = 512;
    case 18, n = 2000; p = 512;
        
    otherwise, return
end

v = 0.5;
k = 5; 
K = 7;
lambda = 0.5*sqrt(1/n)*log(n*p); % MODIFY
SNR = 5;

maxit = 20; tol = 10^-6;

nrun = 30; % MODIFY

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
save(['mat/Stmp_' num2str(version) '.mat'],'J','Ln');
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
epsil = 1e-6;

success_supp = 1e5;

for version = 1:num_versions
    load(['mat/S_' num2str(version) '.mat']);
    nrun = size(Ln,2); suc = 0; 
    for run = 1:nrun
        cur_supp = sum(Ln(:,run) > 1e-5);
        if min(Ln(J(:,run),run)) > epsil && ...
            cur_supp < success_supp
            suc = suc + 1;
        end
        supp(run,version) = cur_supp;
    end
    prob(version) = suc/nrun;
end


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

%% better box plot

samples = [200, 500, 800, 1100, 1500, 2000];
%% box plots
figure;
set(gca,'FontSize',16); 

G = []
v1 = []
v2 = []

% 3 lines
for ii = 1:3
    
    % 6 different sample sizes to try
    for jj = 1:6
        G = [G, supp(:, jj + (6*(ii-1)))']
        % 15 nruns
        v1 = [v1, repmat(samples(jj), 1, 30)]
    end
end

v2 = [repmat({'p=128'}, 1, 180), repmat({'p=256'}, 1, 180), ...
      repmat({'p=512'}, 1, 180)]
  
boxplot(G', {v2'; v1'}, 'factorseparator',1, 'factorgap',5,...
    'colorgroup',v2', 'labelverbosity','majorminor')

title('Number of selected variables')




set(gca, 'units', 'centimeters')
outpos = get(gca, 'OuterPosition')

outpos(4) = outpos(4)*0.8
outpos(3) = outpos(3)*1.4 %width
set(gca, 'OuterPosition', outpos)


set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [outpos(3), outpos(4)])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition', [0,0, outpos(3), outpos(4)])

h = gcf;
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/S_support_box.pdf');



