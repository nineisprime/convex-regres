function [beta,h,obj] = Piecewise_Convex(X,y,lambda,mu,maxit)
[p,n] = size(X); beta = zeros(p,n); S = zeros(n,n); W = zeros(n,n); M = zeros(p,n);
iter = 0; obj = zeros(1,maxit); 
while iter < maxit
    iter = iter + 1;
    
    C = beta-1/mu*M; C = repmat(max(1-lambda/mu./sqrt(sum(C.^2,2)),0),[1,n]).*C;
    
    xb = sum(X.*beta,1);
    h = 1/mu*y-1/mu*(sum(W,2)-sum(W,1)')+sum(S,2)-sum(S,1)'+ ...
        X'*sum(beta,2) + n*xb' - sum(xb) - beta'*sum(X,2);
    h = 1/(1/mu+2*n)*(h+2*mu*sum(h));

    for i = 1:n
        Xi = X - repmat(X(:,i),[1,n]); % Could be computed in advance.
        beta(:,i) = (eye(p)+Xi*Xi')\(C(:,i)+1./mu*M(:,i)+Xi*(h-h(i)-S(:,i)+1/mu*W(:,i)));
        xbi = Xi'*beta(:,i); S(:,i) = max(h-h(i)-xbi+1/mu*W(:,i),0);
        W(:,i) = W(:,i) + mu*(h-h(i)-xbi-S(:,i)); M(:,i) = M(:,i) + mu*(C(:,i)-beta(:,i));
    end
    
    obj(iter) = 0.5*sum((y-h).^2,1) + lambda*sum(sqrt(sum(C.^2,2)));
    mu = min(1.1*mu,10^6);
    if mod(iter,10)==0, 
        disp(['Iter ' num2str(iter) ' Obj ' num2str(obj(iter))]); 
        figure(999); clf; subplot(2,1,1); plot(1:n,y,'ko',1:n,h,'r+');
        subplot(2,1,2); plot(sqrt(sum(beta.^2,2)),'r.'); drawnow;
    end
end
return

%% Testing quadratic function.
clc; clear all; close all; randn('state',0); randn('state',0);

p = 100; n = 200; s = 5; ord = randperm(p); b = zeros(1,p); b(ord(1:s)) = 1;
Xs = randn(n,p)/sqrt(p); Ys = sum(Xs(:,b==1).^2,2);

lambda = 0.001; mu = 0.0001; maxiter = 200;
beta1 = Piecewise_Convex(Xs',Ys, lambda, mu, maxiter); 
beta2 = Piecewise_Convex(Xs',Ys, 0, mu, maxiter);
  
figure(1); subplot(3,1,1); plot(b,'r.'); title('True sparsity pattern');
subplot(3,1,2); plot(sqrt(sum(beta1.^2,2)),'k.'); title('With joint sparse penalty');
subplot(3,1,3); plot(sqrt(sum(beta2.^2,2)),'b.'); title('Without joint sparse penalty');