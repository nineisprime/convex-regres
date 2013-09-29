function [beta,h,obj] = Piecewise_Convex(X,y,lambda,mu,maxit)
[p,n] = size(X); beta = zeros(p,n); S = zeros(n,n); W = zeros(n,n); 
M = zeros(p,n); iter = 0; obj = []; 
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
    
    if mod(iter,50)==0, 
        disp(['   Iter ' num2str(iter) ' Obj ' num2str(obj(iter))]); 
        figure(998); subplot(2,1,1); plot(1:n,y,'ko',1:n,h,'r+');
        subplot(2,1,2); plot(sqrt(sum(beta.^2,2)),'r.'); drawnow;
    end
end
return