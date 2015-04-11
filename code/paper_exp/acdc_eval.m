% IN:
%  "Xt" -- (p--by--nt) test data
%  "pattern" -- (p) either 1 (convex) or -1 (concave)

% OUT:
%   "ht" -- (p--by--nt) predictive components for each test data
%   used by "acdc_boston.m"

% ASSUMPTION:
% all components are non-zero

% FORMULA:
% for a single x_test
% f_test = max_i  h_ik + beta_ik*(x_testk - x_ik)

function ht = acdc_eval(Xt,X,beta,h,pattern)
[p,n] = size(X); 
nt = size(Xt,2); 
ht = zeros(p,nt);

for d = 1:p
    [vv,ss] = sort(X(d,:)); 
    
    if pattern(d) == 1
        ht(d,:) = max(repmat(h(d,ss(1:n-1)),[nt,1]) + repmat(beta(d,:),[nt,1]) ...
            .*(repmat(Xt(d,:)',[1,n-1]) - repmat(X(d,ss(1:n-1)),[nt,1])),[],2);
    else
        ht(d,:) = min(repmat(h(d,ss(1:n-1)),[nt,1]) + repmat(beta(d,:),[nt,1]) ...
            .*(repmat(Xt(d,:)',[1,n-1]) - repmat(X(d,ss(1:n-1)),[nt,1])),[],2);
    end
end


return

