function [beta,alpha,obj,h,c,Ln] = Piecewise_KPlanesADMM(X,y,c,lambda,mu,maxit,beta,alpha)
[p,n] = size(X); K = max(c); 

if (~exist('beta'))
    beta = zeros(p+1,K);  
else
    beta = [alpha; beta];
    beta = beta(:,1:K);
end

M = zeros(p,K); C = zeros(p,K); iter = 0; obj = []; h = zeros(n,1);
while iter < maxit 
    iter = iter + 1;
    
    % C = beta(2:p+1,:)-1/mu*M; C = repmat(max(1-lambda/mu./sqrt(sum(C.^2,2)),0),[1,K]).*C; 
    
    for k = 1:K
        Xk = [ones(1,sum(c==k)); X(:,c==k)]; 
        beta(:,k) = (1/mu*Xk*Xk'+diag([10^-3;ones(p,1)]))\...
                    (1/mu*Xk*y(c==k)+[0;C(:,k)+1/mu*M(:,k)]);
        M(:,k) = M(:,k) + mu*(C(:,k)-beta(2:p+1,k)); h(c==k) = Xk'*beta(:,k);
    end
    
    for j = 1:p, C(j,:) = Piecewise_Inf(beta(j+1,:)-1/mu*M(j,:),lambda/mu); end
    
    Ln = max(abs(C),[],2); obj(iter) = 0.5*sum((y-h).^2,1) + lambda*sum(Ln);
    
    if (iter > 1 && abs(obj(iter) - obj(iter-1)) < obj(iter)*1e-5)
        break;
    end
    
end
alpha = beta(1,:); beta = C; resp = X'*beta + ones(n,1)*alpha; 
[nul,c] = max(resp,[],2); [c1,c2,c] = unique(c);
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