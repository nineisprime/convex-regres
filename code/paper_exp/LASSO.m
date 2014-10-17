function LASSO(data)
clc; close all; format long; randn('state',0); rand('state',0);

switch data
    case 1, load housing.data; 
        X = housing(:,1:end-1)'; 
        y = housing(:,end);
        [p,n] = size(X); 
        q = p; 
    case 2, load diabetes; X = diabetes.x'; y = diabetes.y; [p,n] = size(X); q = p;
    otherwise, return
end

mrun = 10; 
b = 5; 
m = fix(n/b); 
err = zeros(2,q+1,mrun);

for run = 1:mrun
    ord = randperm(n); 
    for i = 1:b
        ndx = ones(n,1)==0; 
        if i<b, 
            ndx(ord((i-1)*m+1:i*m))=1==1; 
        else
            ndx(ord((i-1)*m+1:end))=1==1; 
        end
        erri = LASSO_Main(X,y,ndx,q); 
        
        err(:,:,run) = err(:,:,run) + erri./([b-1; 1]*ones(1,q+1));
    end
end

out = [(0:q)' mean(err,3)' std(err,1,3)']
save(['acdc_boston/lasso_' num2str(data) '.mat'],'out');
return

function [err,yo] = LASSO_Main(X,y,ndx,q)
[p,n] = size(X); 

[X1,X2] = cv_standardize(X,ndx);

y1 = y(~ndx); 
ym = mean(y1); 
y2 = y(ndx);
beta0 = lars(X1',y1-ym,'lasso',0,0,[],0); 
num = sum(abs(beta0)>10^-8,2); 
err = zeros(2,q+1); 
yo = zeros(n,q+1);
for j = 1:q+1 % Refitting
    row = find(num==j-1,1); 
    active = find(abs(beta0(row,:))>10^-8); 
    beta1 = (X1(active,:)')\(y1-ym); % beta1 = beta0(row,active)';
    yo(~ndx,j) = X1(active,:)'*beta1+ym; 
    yo(ndx,j) = X2(active,:)'*beta1+ym;
    err(:,j) = [sum((yo(~ndx,j)-y1).^2); sum((yo(ndx,j)-y2).^2)]/n;
end
return