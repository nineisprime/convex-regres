function [beta,D,h,obj] = Quadratic_Infinity(X,y,lambda,gamma,mu,maxit)
[p,n] = size(X); beta = zeros(p,n); S = zeros(n,n); W = zeros(n,n); 
M = zeros(p,2*n); iter = 0; obj = []; C = zeros(p,2*n); D = zeros(p,n);

num_lambda = length(lambda);
if (num_lambda < p)
    lambda = ones(p,1)*lambda;
end

Delta = zeros(p,n,n);
for j=1:n
    Delta(:,j,:) = (X(:,j)*ones(1,n) - X).^2;
end

while iter < maxit 
    iter = iter + 1;
    for j = 1:p, C(j,:) = Piecewise_Inf([beta(j,:),D(j,:)]-1/mu*M(j,:),lambda(j)/mu); end
    
    Dexpand = permute(repmat(D,[1,1,n]),[1,3,2]);
    tmp = reshape(Delta.*Dexpand,[p,n*n]);
    DDelta = reshape(ones(1,p)*tmp,n,n);
    DDelta = shiftdim(DDelta);
    
    % DDelta(j,i) is D_i'Delta_{ji}, (n--by--n) matrix

    xb = sum(X.*beta,1); %a row of dimension n
    h = 1/mu*y-1/mu*(sum(W,2)-sum(W,1)')+sum(S,2)-sum(S,1)'+ ...
        X'*sum(beta,2) + n*xb' - sum(xb) - beta'*sum(X,2) + ...
        sum(DDelta,2) - sum(DDelta,1)';
    h = 1/(1/mu+2*n)*(h+2*mu*sum(h));
    
    
    %W here is (ji) j indexes row, i indexes column
    
    % X'*sum(beta,2) is a vector, dimension n
    % sum(beta,2) has dimension p
   
    % X'*sum(beta,2)   X'*beta   is n*n matrix, row by X, col by beta
    % sum over columns of X'*beta, produce a column
    
    for i = 1:n
        Xi = X - repmat(X(:,i),[1,n]); 
        beta(:,i) = (eye(p)+Xi*Xi')\(C(:,i)+1./mu*M(:,i)+Xi*(h-h(i)-DDelta(:,i)-S(:,i)+1/mu*W(:,i)));
        xbi = Xi'*beta(:,i); S(:,i) = max(h-h(i)-xbi-DDelta(:,i)+1/mu*W(:,i)-gamma,0);
        %xbi is (n--by--1)
        
        %update for D
        D(:,i) = (eye(p) + Delta(:,:,i)*Delta(:,:,i)')\...
            (C(:,i+n)+1./mu*M(:,i+n)+Delta(:,:,i)*(h-h(i)-xbi-S(:,i)+1/mu*W(:,i)));
        
        W(:,i) = W(:,i) + mu*(h-h(i)-xbi-DDelta(:,i)-S(:,i)); 
        M(:,i) = M(:,i) + mu*(C(:,i)-beta(:,i));
        M(:,i+n) = M(:,i+n) + mu*(C(:,i+n) - D(:,i));
    end
    
    obj(iter) = 0.5*sum((y-h).^2,1) + lambda'*max(abs([beta,D]),[],2) + gamma*sum(sum(S));
    
    if mod(iter,100)==0, 
        disp(['   Iter ' num2str(iter) ' Obj ' num2str(obj(iter)) ' sparsity: ' num2str(sum(max(abs(C),[],2) < 1e-6))]);
        disp(['admm error: ' num2str(norm(C-[beta,D],'fro'))]);
        figure(999); subplot(2,1,1); plot(1:n,y,'ko',1:n,h,'r+');
        subplot(2,1,2); plot(max(abs(C),[],2),'r.'); drawnow;
    end
        
end

beta = C(:,1:n);
D = C(:,(n+1):(2*n));

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