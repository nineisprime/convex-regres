
% ** DEPRECATED **
%  
% 
% 

% IN:
%  "Xt" -- (p--by--nt) test data

% OUT:
%   "ht" -- (p--by--nt) predictive components for each test data

% Formula:
% for a single x_test
% f_test = max_i  h_ik + beta_ik*(x_testk - x_ik)

function ht = SCAM_Eval(Xt,X,beta,h)
[p,n] = size(X); 
nt = size(Xt,2); 
ht = zeros(p,nt);

for d = 1:p
    [vv,ss] = sort(X(d,:)); 
    ht(d,:) = max(repmat(h(d,ss(1:n-1)),[nt,1]) + repmat(beta(d,:),[nt,1]) ...
        .*(repmat(Xt(d,:)',[1,n-1]) - repmat(X(d,ss(1:n-1)),[nt,1])),[],2);
end


return

