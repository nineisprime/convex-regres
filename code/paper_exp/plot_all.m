


%% Chapter 1.
% plot predictive error / errorbar (y-axis) against
% number of selected variables (x-axis)

% depends on "acdc_boston/...mat" files produced by
% "acdc_boston.m" and
% "LASSO_boston.m"

data = 1;
clc; close all; format long; randn('state',0); rand('state',0);

lambda0 = {[0.00 .005 0.01 .015 0.02 0.04 0.06 ...
            0.07 0.08 0.09 0.10 0.12 0.20 0.30 0.40 ...
            0.50, 0.60, 0.70, 0.80, 0.90, 1, 1.2, 1.4 ...
            1.6, 1.8 4],...
           };
       
nd = size(lambda0,2); 
nv = zeros(1,nd); 

for i = 1:nd, nv(i) = length(lambda0{i}); end


switch data
    case 1, load housing.data; 
        X = housing(:,1:end-1)'; 
        y = housing(:,end);
    % case 2, 
    otherwise, return
end

[p,n] = size(X); 
acdc_out = zeros(nv(data),5);

for version = sum(nv(1:data-1))+1:sum(nv(1:data))
    if (version == 1)
        continue
    end
    load(['acdc_boston/acdc_' num2str(version) '.mat']); % loads "out" var 
    acdc_out(version-sum(nv(1:data-1)),:) = out;
end


output(1,:) = [];
%output(end-1,:) = []
%output(end-1,:) = []
%output(end,1) = 0

load(['acdc_boston/lasso_' num2str(data) '.mat']);
figure(1); set(gca,'FontSize',18);
lin_out = out;

errorbar(lin_out(:,1), lin_out(:,3), lin_out(:,5),'g--','LineWidth',2); hold on

%[~,ord] = sort(output(:,1), 'descend');
errorbar(acdc_out(:,1), acdc_out(:,3), acdc_out(:,5),'k:','LineWidth',2); hold off

legend('LASSO','ACDC'); 
xlim([-1 out(end,1)+1])
xlabel('Average number of seleted features'); 
ylabel('Testing Mean Squared Error')
title('Testing errors for LASSO and ACDC');


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
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/MSEacdc.pdf');



%% Chapter 2, plot regularization path
% plot acdc regularization path
% depends on "acdc_QP"


% Analysis using all the data.
names = {{'CRIM','ZN','INDUS','CHAS','NOX','RM','AGE','DIS','RAD','TAX','PTRATIO','B','LSTAT'},...
         }; 

X1 = cv_standardize(X,~ones(1,n)); 
ym = mean(y); 
maxit = 20; 
tol = 10^-6; 
beta = []; 
Ln = [];
for run = 1:nv(data)
    [beta(:,:,run),h,obj,Lnvex,Lncave] = acdc_QP(X1,y-ym,lambda0{data}(run),maxit,tol); 
    Ln(:,run) = max(Lnvex, Lncave);
end
ratio = sum(Ln,1)/sum(Ln(:,1),1);

figure(2); 
set(gca,'FontSize',18); 
set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0 0 0;1 1 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.'); 
fg2 = plot(ratio,Ln,'LineWidth',2);
legend(names{data}); 
xlabel('Normalized |f_j|_{\infty,1}'); 
ylabel('|f_j|_\infty'); 
ylim([0,27])
title('ACDC')


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
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/acdc_path.pdf');


% Plot lasso regularization path
% dependent on "lars.m" only

betao = lars(X1',y-ym,'lasso',0,0,[],0); 
ratioo = sum(abs(betao),2)/sum(abs(betao(end,:)));
figure(3); 
set(gca,'FontSize',18); 
fg3 = plot(ratioo,abs(betao),'LineWidth',2); 
legend(names{data}); 
xlabel('Normalized |\beta|_{1}'); 
ylabel('\beta'); 
title('LASSO')


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
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/lasso_path.pdf');



%% Chapter 3. plot example additive component 

%figure(4); 
%set(gca,'FontSize',12); 
%fg4 = plot(ratio,squeeze(mean(beta,2)),'LineWidth',2); 
%legend(names{data}); 
%xlabel('Normalized |\beta|_{\infty,1}'); 
%ylabel('mean(\beta)'); 
%title('ACDC')

% A special choice of lambda.
%[beta0,h0,obj0,Ln0] = SCAM_QP(X1,y-ym,lambda0{data}(10),maxit,tol); 
my_lambda = 0.7;
[beta0,h0,obj0,Lnvex0,Lncave0] = acdc_QP(X1,y-ym, my_lambda ,maxit,tol); 
Ln0 = max(Lnvex0, Lncave0);
active = find(abs(Ln0)>10^-6); 
num = length(active);

pattern( abs(Lnvex0)>1e-8 & Lnvex0 >= Lncave0) = 1;
pattern( abs(Lncave0)>1e-8 & Lnvex0 < Lncave0) = -1;
% eliminate 0
pattern = pattern(active);
pattern
[beta1,h1,obj1] = acdc_refit(X1(active,:),y-ym,pattern,maxit,tol); 

len = 1000; xd = zeros(num,len);

for d = 1:num 
    xd(d,:) = linspace(min(X1(active(d),:)),max(X1(active(d),:)),len); 
end

hd = acdc_eval(xd,X1(active,:),beta1,h1,pattern); figure(5);

for d = 1:num

    subplot(2,2,d); 
    set(gca,'FontSize',14);  % ,X1(active(d),:),h1(d,:),'k.'
    plot(xd(d,:),hd(d,:),'r-','LineWidth',2); 
    xlim([xd(d,1) xd(d,end)]); 
    xlabel('Feature'); 
    ylabel('Response'); 
    ylim([-14,20])
    title(names{data}{active(d)});
end

ti = get(gca,'TightInset');
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

h = gcf;
saveas(h, '/Users/minxu/dropbox/minx/research/convex_regr/tex/scam/figs/acdc_functs.pdf');



return