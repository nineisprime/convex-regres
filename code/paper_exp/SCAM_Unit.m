
% ** DEPRECATED **
% replaced by "cv_standardize.m"

% OUT:
% "X1 = transformed X(:,ndx) "
% "X2 = transformed X(:,~ndx)"

% standardizes X1, transform X2 accordingly

% a predictor trained on standardized X1 can
% be used on the returned X2

function [X1,X2,Xm,Xstd] = SCAM_Unit(X,ndx)

[p,n] = size(X); 
n1 = sum(~ndx); 

X1 = X(:,~ndx); 
Xm = mean(X1,2);

X1 = X1 - Xm*ones(1,n1); 
Xstd = sqrt(sum(X1.^2,2))+eps; 
X1 = X1./(Xstd*ones(1,n1)); 

X2 = (X(:,ndx) - Xm*ones(1,n-n1))./(Xstd*ones(1,n-n1));
return
