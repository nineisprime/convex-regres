function R(version)
% Hessian changes
% Non-diagonal quadratic with noise.
% Independent


clc; close all; format long; randn('state',0); rand('state',0); 
switch version
    case 00, n = 100; alpha = 0.2;
    case 01, n = 100; alpha = 0.0;
    case 02, n = 200; alpha = 0.0;
    case 03, n = 300; alpha = 0.0;
    case 04, n = 400; alpha = 0.0;
    case 05, n = 500; alpha = 0.0;
    case 06, n = 600; alpha = 0.0;
    case 07, n = 700; alpha = 0.0;
    case 08, n = 800; alpha = 0.0;
    case 09, n = 900; alpha = 0.0;
    case 10, n =1000; alpha = 0.0;
        
    case 11, n = 100; alpha = 0.2;
    case 12, n = 200; alpha = 0.2;
    case 13, n = 300; alpha = 0.2;
    case 14, n = 400; alpha = 0.2;
    case 15, n = 500; alpha = 0.2;
    case 16, n = 600; alpha = 0.2;
    case 17, n = 700; alpha = 0.2;
    case 18, n = 800; alpha = 0.2;
    case 19, n = 900; alpha = 0.2;
    case 20, n =1000; alpha = 0.2;
        
    case 21, n = 100; alpha = 0.5;
    case 22, n = 200; alpha = 0.5;
    case 23, n = 300; alpha = 0.5;
    case 24, n = 400; alpha = 0.5;
    case 25, n = 500; alpha = 0.5;
    case 26, n = 600; alpha = 0.5;
    case 27, n = 700; alpha = 0.5;
    case 28, n = 800; alpha = 0.5;
    case 29, n = 900; alpha = 0.5;
    case 30, n =1000; alpha = 0.5;
        
    case 31, n = 100; alpha = 1.0;
    case 32, n = 200; alpha = 1.0;
    case 33, n = 300; alpha = 1.0;
    case 34, n = 400; alpha = 1.0;
    case 35, n = 500; alpha = 1.0;
    case 36, n = 600; alpha = 1.0;
    case 37, n = 700; alpha = 1.0;
    case 38, n = 800; alpha = 1.0;
    case 39, n = 900; alpha = 1.0;
    case 40, n =1000; alpha = 1.0;
    otherwise, return
end

p = 128; k = 5; 
lambda = 0.5*sqrt(1/n)*log(n*p); % MODIFY
maxit = 20; tol = 10^-6;
nrun = 100; 
J = ones(p,nrun)==0; Ln = zeros(p,nrun);

if alpha<10^-16
    Q = eye(k); 
else
    Q = 0.5*eye(k); 
    B0 = rand(k,k) < alpha;
    for i = 1:k
        for j = 1:i
            if B0(i,j) && i~=j
                Q(i,j) = 0.5; 
            end; 
        end; 
    end; Q = Q + Q';
end

for run = 1:nrun
    ord = randperm(p); J(ord(1:k),run) = 1==1; 
    X = randn(p,n); 
    y = sum(X(J(:,run),:).*(Q*X(J(:,run),:)),1)' + randn(n,1);
    % main function call
    [beta,h,obj,Lnvex,Lncave] = acdc_QP(X, y-mean(y), lambda, maxit, tol);
    
    Ln(:,run) = max(Lnvex,Lncave);
    %[beta,h,obj,Ln(:,run)] = SCAM_QP(X,y-mean(y),lambda,maxit,tol);
    disp(['version=' num2str(version) ' run=' num2str(run)]);
end
save(['R_' num2str(version) '.mat'],'J','Ln'); %MODIFY
%save('tmp.mat', 'J', 'Ln');
return

%% Reading the result:
clc; clear all; close all; 
alpha = [0.0 0.2 0.5 1.0]; k = 5;
for v = 1:length(alpha)
    randn('state',0); rand('state',0);
    if alpha(v)<10^-16
        Q = eye(k); 
    else
        Q = 0.5*eye(k); 
        B0 = rand(k,k) < alpha(v);
        
        for i = 1:k 
            for j = 1:i 
                if B0(i,j) && i~=j
                    Q(i,j) = 0.5; 
                end; 
            end; 
        end; 
    
        Q = Q + Q'; 
    end
    figure(1); subplot(2,2,v); imagesc(Q); 
    title(['\alpha=' num2str(alpha(v))]); 
    caxis([0 1.5]);colorbar;
end


prob = zeros(1,40); epsil = 10^-6;
for version = 1:40
    load(['R_' num2str(version) '.mat']);
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
legend('\alpha=0.0','\alpha=0.2','\alpha=0.5','\alpha=1.0');
xlabel('Number of samples'); ylabel('Probability of recovery');
title('Probability of support recovery');


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
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/Curve2.pdf')