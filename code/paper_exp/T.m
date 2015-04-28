function T(version)
% p changes
% function Diagonal quadratic with noise.
% independent gaussian

clc; close all; format long; randn('state',0); rand('state',0); 
switch version
    case 00, n = 100; p = 064;
    case 01, n = 100; p = 064;
    case 02, n = 200; p = 064;
    case 03, n = 300; p = 064;
    case 04, n = 400; p = 064;
    case 05, n = 500; p = 064;
    case 06, n = 600; p = 064;
    case 07, n = 700; p = 064;
    case 08, n = 800; p = 064;
    case 09, n = 900; p = 064;
    case 10, n =1000; p = 064;
        
    case 11, n = 100; p = 128;
    case 12, n = 200; p = 128;
    case 13, n = 300; p = 128;
    case 14, n = 400; p = 128;
    case 15, n = 500; p = 128;
    case 16, n = 600; p = 128;
    case 17, n = 700; p = 128;
    case 18, n = 800; p = 128;
    case 19, n = 900; p = 128;
    case 20, n =1000; p = 128;
        
    case 21, n = 100; p = 256;
    case 22, n = 200; p = 256;
    case 23, n = 300; p = 256;
    case 24, n = 400; p = 256;
    case 25, n = 500; p = 256;
    case 26, n = 600; p = 256;
    case 27, n = 700; p = 256;
    case 28, n = 800; p = 256;
    case 29, n = 900; p = 256;
    case 30, n =1000; p = 256;
  
    case 31, n = 100; p = 512;
    case 32, n = 200; p = 512;
    case 33, n = 300; p = 512;
    case 34, n = 400; p = 512;
    case 35, n = 500; p = 512;
    case 36, n = 600; p = 512;
    case 37, n = 700; p = 512;
    case 38, n = 800; p = 512;
    case 39, n = 900; p = 512;
    case 40, n =1000; p = 512;
    otherwise, return
end

k = 5; 
lambda = 0.5*sqrt(1/n)*log(n*p); %MODIFY

alpha = 0.5;

maxit = 20; tol = 10^-6;
nrun = 100; J = ones(p,nrun)==0; 
Ln = zeros(p,nrun);



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
    
    ord = randperm(p); J(ord(1:k),run) = 1==1; X = randn(p,n); 
    
    y = sum(X(J(:,run),:).*(Q*X(J(:,run),:)),1)' + randn(n,1);
    
    [beta,h,obj,Lnvex,Lncave] = acdc_QP(X,y-mean(y),lambda,maxit,tol);
    Ln(:,run) = max(Lnvex, Lncave);
    disp(['version=' num2str(version) ' run=' num2str(run)]);
end
save(['mat/T_' num2str(version) '.mat'],'J','Ln');
return

%% Reading the result:
clc; clear all; close all; 
prob = zeros(1,40); epsil = 10^-6;
for version = 1:40
    load(['mat/T_' num2str(version) '.mat']);
    nrun = size(Ln,2); suc = 0; 
    for run = 1:nrun
        if max(Ln(~J(:,run),run)) < epsil && min(Ln(J(:,run),run)) > epsil
            suc = suc + 1;
        end
    end
    prob(version) = suc/nrun;
end

figure(2); set(gca,'FontSize',14); 
plot(100:100:1000,prob(1:10),'r.-',100:100:1000,prob(11:20),'b.-',...
    100:100:1000,prob(21:30),'g.-',100:100:1000,prob(31:40),'k.-','LineWidth',2);
legend('p=64','p=128','p=256','p=512','Location','SouthEast');
xlabel('Number of Samples'); 
ylabel('Probability of Recovery');
title('Probability of Recovery');


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
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/CurveT.pdf')


