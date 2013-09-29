function [beta,z,obj,Ln] = Piecewise_Additive(X,y,lambda,mu,maxit)
[p,n] = size(X); beta = zeros(p,n); z = zeros(p,n);
S = zeros(p,n,n); W = zeros(p,n,n); M = zeros(p,n); C = zeros(p,n); q = zeros(p,1);
iter = 0; obj = []; tmp0 = 1/mu*repmat(y',[p,1]);
while iter < maxit
    iter = iter + 1;
    for j = 1:p, C(j,:) = Piecewise_Inf(beta(j,:)-1/mu*M(j,:),lambda/mu); end
    
    tmp1 = repmat(sum(beta,2),[1,n]).*X - repmat(sum(beta.*X,2),[1,n]) ...
        - repmat(sum(X,2),[1,n]).*beta + n*X.*beta + squeeze(sum(S,3))-squeeze(sum(S,2)) ...
        - 1/mu*(squeeze(sum(W,3))-squeeze(sum(W,2))) - 1/mu*repmat(q,[1,n]);
    z = tmp0 + tmp1;
    [ialpha,igamma,is] = Piecewise_Inverse(2*n*ones(p,1),ones(p,1),-1/mu);
    z = repmat(ialpha,[1,n]).*z - is*igamma*(igamma'*z); zs = sum(z,2);
    [alpha,gamma,s] = Piecewise_Inverse(ones(p,1)-n*ialpha,igamma,-n*is);
    zs = alpha.*zs - s*gamma*(gamma'*zs); zs = ialpha.*zs - is*igamma*(igamma'*zs);
    z = z + repmat(zs,[1 n]);
    
    for i = 1:n
        Xi = X - repmat(X(:,i),[1,n]); zi = z - repmat(z(:,i),[1,n]);
        beta(:,i) = (C(:,i)+1/mu*M(:,i)+sum(Xi.*(zi-S(:,:,i)+1/mu*W(:,:,i)),2))./(1+sum(Xi.^2,2));
        xbi = Xi.*repmat(beta(:,i),[1,n]);
        S(:,:,i) = max(zi-xbi+1/mu*W(:,:,i),0); W(:,:,i) = W(:,:,i) + mu*(zi-xbi-S(:,:,i));
    end
    q = q + mu*sum(z,2); M = M + mu*(C-beta);
    Ln = max(abs(C),[],2); obj(iter) = 0.5*sum((y'-sum(z,1)).^2) + lambda*sum(Ln);
    if mod(iter,10)==0,
        disp(['   Iter ' num2str(iter) ' Obj ' num2str(obj(iter))]);
        figure(1); subplot(2,1,1); plot(1:n,y','ko',1:n,sum(z,1),'r+'); title('Response');
        subplot(2,1,2); plot(Ln,'r.'); title('L\infty norm'); drawnow;
    end
end
return

%% Testing case:
clc; clear all; close all; randn('state',1); rand('state',0);
p = 40; n = 400; k = 3; 
X = randn(p,n)/4; 
%ord = randperm(p); 
ord = 1:p;

y = sum(X(ord(1:k),:).^2,1)';  
A = randn(k,k);
A = A*A'+ 0.5*eye(k);
%A = A/norm(A,2);
A = eye(k);
y = diag(X(ord(1:k),:)'*A*X(ord(1:k),:));

lambda = n*0.01; mu = 0.001; maxit = 1000;
[beta,z,obj,Ln] = Piecewise_Additive(X,y-mean(y),lambda,mu,maxit);