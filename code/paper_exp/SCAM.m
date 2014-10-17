% Run SCAM on Boston housing data. 
% Use 5-fold CV and plot prediction error

function SCAM(version)
clc; close all; format long; randn('state',0); rand('state',0);
lambda0 = {[0.00 .005 0.01 .015 0.02 0.04 0.06 ...
            0.07 0.08 0.09 0.10 0.12 0.20 0.30 0.40],...
           };
nd = size(lambda0,2); 
nv = zeros(1,nd); 

for i = 1:nd, nv(i) = length(lambda0{i}); end

data = sum(version>cumsum(nv)) + 1; 
lambda = lambda0{data}(version - sum(nv(1:data-1)));

switch data
    case 1, load housing.data; X = housing(:,1:end-1)'; y = housing(:,end);
    % case 2, 
    otherwise, return
end

[p,n] = size(X); 
mrun = 10; 
b = 5; 
m = fix(n/b); % 5-fold CV

err = zeros(2,mrun); 
num = zeros(1,mrun);

for run = 1:mrun
    ord = randperm(n);
    for i = 1:b
        ndx = zeros(n,1);
        if i<b
            ndx(ord((i-1)*m+1:i*m))= 1; 
        else
            ndx(ord((i-1)*m+1:end))= 1; 
        end
        
        [erri,numi] = SCAM_Main(X,y,ndx,lambda);
        
        err(:,run) = err(:,run) + erri./[b-1; 1]; 
        num(run) = num(run) + numi/b;
        
        disp(['run=' num2str(run) ' cv=' num2str(i)]);
    end
end

out = [mean(num) mean(err,2)' std(err,1,2)']
save(['SCAM/SCAM_' num2str(version) '.mat'],'out');
return

function [err,num,yo] = SCAM_Main(X,y,ndx,lambda)
[p,n] = size(X); 
[X1,X2] = SCAM_Unit(X,ndx);

y1 = y(~ndx); 
ym = mean(y1); 
y2 = y(ndx); 

maxit = 20; tol = 10^-6;

[beta0,h0,obj0,Ln0] = SCAM_QP(X1,y1-ym,lambda,maxit,tol);

active = find(abs(Ln0)>10^-8); 
num = length(active);

[beta1,h1,obj1,Ln1] = SCAM_QP(X1(active,:),y1-ym,0,maxit,tol);

h2 = SCAM_Eval(X2(active,:),X1(active,:),beta1,h1);

yo = zeros(n,1); 
yo(~ndx) = sum(h1,1)'+ym; 
yo(ndx) = sum(h2,1)'+ym;

err = [sum((yo(~ndx)-y1).^2); sum((yo(ndx)-y2).^2)]/n;
return