% Run ACDC on Boston housing data. 
% Use 5-fold CV and save prediction error
%


% acdc analogue of the "SCAM.m" file in same directory

% OUT:
%  a tuple, saved to "acdc_boston/" directory
%  data used by "plot.m"


function acdc_boston(version)

clc; close all; format long; randn('state',0); rand('state',0);

% potentially many sequences in lambda0, each sequence for a particular
% data set.
% version 1,...,m1 correspond to data 1 and lambda for data 1
% version m1+1,....,m1+m2 correspond to data 2 etc.

lambda0 = {[0.00 .005 0.01 .015 0.02 0.04 0.06 ...
            0.07 0.08 0.09 0.10 0.12 0.20 0.30 0.40 ...
            0.50, 0.60, 0.70, 0.80, 0.90, 1, 1.2, 1.4 ...
            1.6, 1.8, 4],...
           };

       
nd = size(lambda0,2); % possibly many different lambda-sequences
nv = zeros(1,nd); 

for i = 1:nd, 
    nv(i) = length(lambda0{i}); 
end

data = sum(version>cumsum(nv)) + 1; 

lambda = lambda0{data}(version - sum(nv(1:data-1)));

switch data
    case 1, load housing.data; 
        X = housing(:,1:end-1)'; 
        y = housing(:,end);
    % case 2, 
    otherwise, return
end

[p,n] = size(X); 
mrun = 3; 
b = 5; 
m = fix(n/b); % 5-fold CV

err = zeros(2,mrun); 
num = zeros(1,mrun);

for run = 1:mrun
    ord = randperm(n);
    
    %iterate over folds
    for i = 1:b
        ndx = zeros(n,1);
        
        % set validation data
        if i<b
            ndx(ord((i-1)*m+1:i*m))= 1; 
        else
            ndx(ord((i-1)*m+1:end))= 1; 
        end
        
       
        [erri,numi] = acdc_main(X,y,ndx,lambda);
        erri
        numi
        
        err(:,run) = err(:,run) + erri./[b-1; 1]; 
        num(run) = num(run) + numi/b;
        
        disp(['run=' num2str(run) ' cv=' num2str(i)]);
    end
end

out = [mean(num) mean(err,2)' std(err,1,2)'];
save(['acdc_boston/acdc_' num2str(version) '.mat'],'out');
return



function [err,num,yo] = acdc_main(X,y,ndx,lambda)
[p,n] = size(X); 
[X1,X2] = cv_standardize(X,ndx);

y1 = y(~ndx); 
ym = mean(y1); 
y2 = y(ndx==1); 

maxit = 20; tol = 10^-6;

[beta0, h0, obj0, Lnvex0, Lncave0] = acdc_QP(X1,y1-ym,lambda,maxit,tol);

active = find(max(abs(Lnvex0), abs(Lncave0))>1e-7); 
num = length(active);

%num
% pattern is 1 if Lnvex0 is non-zero and large
%           -1 if Lncave0 is non-zero and large
% else 0
pattern( abs(Lnvex0)>1e-8 & Lnvex0 >= Lncave0) = 1;
pattern( abs(Lncave0)>1e-8 & Lnvex0 < Lncave0) = -1;
% eliminate 0
pattern = pattern(active);
pattern

%if (num ~= length(pattern))
%    Lnvex0
%    Lncave0
%end

[beta1,h1,obj1] = acdc_refit(X1(active,:),y1-ym,pattern,maxit,tol);

% currently here
h2 = acdc_eval(X2(active,:),X1(active,:),beta1,h1,pattern);

yo = zeros(n,1); 
yo(~ndx) = sum(h1,1)'+ym; 
yo(ndx==1) = sum(h2,1)'+ym;

err = [sum((yo(~ndx)-y1).^2); sum((yo(ndx==1)-y2).^2)]/n;
return