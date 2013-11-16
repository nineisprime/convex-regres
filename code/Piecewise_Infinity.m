function [beta,h,obj] = Piecewise_Infinity(X,y,lambda,mu,maxit)
[p,n] = size(X); beta = zeros(p,n); S = zeros(n,n); W = zeros(n,n); 
M = zeros(p,n); iter = 0; obj = []; C = zeros(p,n);

num_lambda = length(lambda);
if (num_lambda < p)
    lambda = ones(p,1)*lambda;
end

while iter < maxit 
    iter = iter + 1;
    for j = 1:p, C(j,:) = Piecewise_Inf(beta(j,:)-1/mu*M(j,:),lambda(j)/mu); end
    
    xb = sum(X.*beta,1); %a row of dimension n
    h = 1/mu*y-1/mu*(sum(W,2)-sum(W,1)')+sum(S,2)-sum(S,1)'+ ...
        X'*sum(beta,2) + n*xb' - sum(xb) - beta'*sum(X,2);
    h = 1/(1/mu+2*n)*(h+2*mu*sum(h));
    
    prev_beta = beta;
    all_xb = zeros(n,n);
    all_prev_xb = zeros(n,n);
    h_mat = zeros(n,n);
    
    for i = 1:n
        Xi = X - repmat(X(:,i),[1,n]); 
        beta(:,i) = (eye(p)+Xi*Xi')\(C(:,i)+1./mu*M(:,i)+Xi*(h-h(i)-S(:,i)+1/mu*W(:,i)));
        xbi = Xi'*beta(:,i); S(:,i) = max(h-h(i)-xbi+1/mu*W(:,i),0); 
        W(:,i) = W(:,i) + mu*(h-h(i)-xbi-S(:,i)); M(:,i) = M(:,i) + mu*(C(:,i)-beta(:,i)); 
        
        %for stopping criteria
        all_xb(i,:) = xbi';
        all_prev_xb(i,:) = prev_beta(:,i)'*Xi;
        h_mat(:,i) = h-h(i)-xbi;
    end
    
    obj(iter) = 0.5*sum((y-h).^2,1);
    
    
    
    %=- verify stopping criteria -=
    criterion1 = norm(all_xb - all_prev_xb,'fro')/n;
    criterion2 = norm(S - h_mat, 'fro')/n;
    criterion3 = norm(C - beta,'fro')/sqrt(n*p);
    
    threshold = 5e-5;
    stop_score = max([criterion1, criterion2, criterion3]);
    if (stop_score < threshold)
        disp('ADMM converges.');
        break;
    end
    
    if mod(iter,10)==0, 
        [criterion1, criterion2, criterion3]
        disp(['   Iter ' num2str(iter) ' Obj ' num2str(obj(iter)) ' Score ' num2str(stop_score)]); 
        figure(999); subplot(2,1,1); plot(1:n,y,'ko',1:n,h,'r+');
        subplot(2,1,2); plot(max(abs(beta),[],2),'r.'); drawnow;
    end    
end
beta = C;
return

function z = Piecewise_Inf(x,v)
% L-infinity shrinkage.
n = length(x); [xa,ndx] = sort(abs(x),'descend');
b = [xa(2:n) 0] - (cumsum(xa)-v)./(1:n); z = zeros(1,n); 
if b(n) < 0
    pos = find(b<0,1); zinf = (sum(abs(x(ndx(1:pos))))-v)/pos;
    z(ndx(1:pos)) = sign(x(ndx(1:pos))).*zinf;
    z(ndx(pos+1:n)) = x(ndx(pos+1:n));
end
return