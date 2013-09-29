function [beta,h,obj,Ln] = Piecewise_InfinityADMM(X,y,lambda,mu,maxit)
[p,n] = size(X); beta = zeros(p,n); S = zeros(n,n); W = zeros(n,n); 
M = zeros(p,n); iter = 0; obj = []; C = zeros(p,n);
while iter < maxit 
    iter = iter + 1;
    for j = 1:p, C(j,:) = Piecewise_Inf(beta(j,:)-1/mu*M(j,:),lambda/mu); end
    
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
    
    Ln = max(abs(C),[],2); obj(iter) = 0.5*sum((y-h).^2,1) + lambda*sum(Ln);
    if mod(iter,10000)==0, 
        disp(['   Iter ' num2str(iter) ' Obj ' num2str(obj(iter))]); 
        figure(1); subplot(2,1,1); plot(1:n,y,'ko',1:n,h,'r+');
        subplot(2,1,2); plot(Ln,'r.'); drawnow;
    end
end
return

function z = Piecewise_Inf(x,v)
n = length(x); [xa,ndx] = sort(abs(x),'descend');
b = [xa(2:n) 0] - (cumsum(xa)-v)./(1:n); z = zeros(1,n); 
if b(n) < 0
    pos = find(b<0,1); zinf = (sum(abs(x(ndx(1:pos))))-v)/pos;
    z(ndx(1:pos)) = sign(x(ndx(1:pos))).*zinf;
    z(ndx(pos+1:n)) = x(ndx(pos+1:n));
end
return