data = 1;
clc; close all; format long; randn('state',0); rand('state',0);
lambda0 = {[0.00 .005 0.01 .015 0.02 0.04 0.06 0.07 0.08 0.09 0.10 0.12 0.20 0.30 0.40],...
           };
nd = size(lambda0,2); nv = zeros(1,nd); for i = 1:nd, nv(i) = length(lambda0{i}); end
switch data
    case 1, load housing.data; X = housing(:,1:end-1)'; y = housing(:,end);
    % case 2, 
    otherwise, return
end

[p,n] = size(X); output = zeros(nv(data),5);
for version = sum(nv(1:data-1))+1:sum(nv(1:data))
    load(['SCAM/SCAM_' num2str(version) '.mat']); output(version-sum(nv(1:data-1)),:) = out;
end
load(['SCAM/LASSO_' num2str(data) '.mat']);
figure(1); set(gca,'FontSize',14);
errorbar(out(:,1),out(:,3),out(:,5),'g--','LineWidth',2); hold on
errorbar(output(:,1),output(:,3),output(:,5),'k:','LineWidth',2); hold off
legend('LASSO','SCAM'); xlim([-1 out(end,1)+1])
xlabel('Average number of seleted features'); ylabel('Testing Mean Squared Error')
title('Testing errors for LASSO and SCAM');

% Analysis using all the data.
names = {{'CRIM','ZN','INDUS','CHAS','NOX','RM','AGE','DIS','RAD','TAX','PTRATIO','B','LSTAT'},...
         }; 

X1 = SCAM_Unit(X,ones(1,n)==0); ym = mean(y); maxit = 20; tol = 10^-6; beta = []; Ln = [];
for run = 1:nv(data)
    [beta(:,:,run),h,obj,Ln(:,run)] = SCAM_QP(X1,y-ym,lambda0{data}(run),maxit,tol); 
end
ratio = sum(Ln,1)/sum(Ln(:,1),1);
figure(2); set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;0 0 0;1 1 0],...
      'DefaultAxesLineStyleOrder','-|--|:|-.'); set(gca,'FontSize',12);
fg2 = plot(ratio,Ln,'LineWidth',2);
legend(names{data}); xlabel('Normalized |\beta|_{\infty,1}'); ylabel('|\beta|_\infty'); title('SCAM')

betao = lars(X1',y-ym,'lasso',0,0,[],0); ratioo = sum(abs(betao),2)/sum(abs(betao(end,:)));
figure(3); set(gca,'FontSize',12); fg3 = plot(ratioo,betao,'LineWidth',2); 
legend(names{data}); xlabel('Normalized |\beta|_{1}'); ylabel('\beta'); title('LASSO')

figure(4); set(gca,'FontSize',12); fg4 = plot(ratio,squeeze(mean(beta,2)),'LineWidth',2); 
legend(names{data}); xlabel('Normalized |\beta|_{\infty,1}'); ylabel('mean(\beta)'); title('SCAM')

% A special choice of lambda.
[beta0,h0,obj0,Ln0] = SCAM_QP(X1,y-ym,lambda0{data}(10),maxit,tol); 
active = find(abs(Ln0)>10^-8); num = length(active);
[beta1,h1,obj1,Ln1] = SCAM_QP(X1(active,:),y-ym,0,maxit,tol); 

len = 1000; xd = zeros(num,len);
for d = 1:num, xd(d,:) = linspace(min(X1(active(d),:)),max(X1(active(d),:)),len); end
hd = SCAM_Eval(xd,X1(active,:),beta1,h1); figure(5); 
for d = 1:num
    subplot(2,2,d); set(gca,'FontSize',12);  % ,X1(active(d),:),h1(d,:),'k.'
    plot(xd(d,:),hd(d,:),'r-','LineWidth',2); xlim([xd(d,1) xd(d,end)]); 
    xlabel('Feature'); ylabel('Response'); title(names{data}{active(d)});
end
return