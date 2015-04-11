% INPUT:
%   X must be (n--by--p)
%   s is the sparsity level
%   H is a PSD matrix, optional

% OUTPUT:
%   Y is (n--by--1)

function [Y] = quadraticForm(X, s, H)

%[n,p] = size(X);

if (~exist('H', 'var'))
   H = randn(s,s);
   H = H'*H;
   D = diag(diag(H));
   H = D^(-1/2)*H*D^(-1/2);
    
end

Y = diag(X(:,1:s)*H*X(:,1:s)');