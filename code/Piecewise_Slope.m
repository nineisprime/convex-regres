function [C,h,obj,Ln] = Piecewise_Slope(X,y,lambda,mu,maxit)
[p,n] = size(X); dX = zeros(p,n-1); ord = zeros(p,n);
for d = 1:p, [Xd,ord(d,:)] = sort(X(d,:)); dX(d,:) = Xd(2:n)-Xd(1:n-1); end
beta = zeros(p,n-1); M = zeros(p,n-1); S = zeros(p,n-2); W = zeros(p,n-2);
C = zeros(p,n-1); O = zeros(p,n-1); q = zeros(p,1); iter = 0; obj = [];
G = gallery('tridiag',-ones(1,n-2),[2 3*ones(1,n-3) 2],-ones(1,n-2));
A = gallery('tridiag',-ones(1,n-1),[1 2*ones(1,n-2) 1],-ones(1,n-1));
B = inv(A+ones(n,n)); B1 = zeros(n,n); B2 = eye(n);
for d = 1:p, B1(ord(d,:),:) = B; B1(:,ord(d,:)) = B1; B2 = B2+1/mu*B1; end

while iter < maxit
    iter = iter + 1;
    for d = 1:p, C(d,:) = Piecewise_Inf(beta(d,:)-1/mu*O(d,:),lambda/mu); end
    
    h = zeros(p,n); eta = 1/mu*repmat(y',[p 1])-1/mu*repmat(q,[1 n]);
    for d = 1:p
        resp = dX(d,:).*beta(d,:)-1/mu*M(d,:);
        eta(d,ord(d,1:n)) = eta(d,ord(d,1:n)) + ([0 resp] - [resp 0]);
        h(d,ord(d,:)) = eta(d,ord(d,:))*B;
    end
    vec = B2\(sum(h,1)');
    for d = 1:p, h(d,ord(d,:)) = h(d,ord(d,:)) - 1/mu*vec(ord(d,:))'*B; end
    
    for d = 1:p
        dh = h(d,ord(d,2:n))-h(d,ord(d,1:n-1));
        gd = C(d,:)+1/mu*O(d,:) + [0 S(d,:)-1/mu*W(d,:)] - [S(d,:)-1/mu*W(d,:) 0] + ...
            dX(d,:).*(dh+1/mu*M(d,:));
        Gd = G; rg = (0:n-2)*(n-1)+(1:n-1); Gd(rg) = Gd(rg) + dX(d,:).^2;
        beta(d,:) = Gd\(gd'); dbeta = beta(d,2:n-1)-beta(d,1:n-2);
        S(d,:) = max(dbeta+1/mu*W(d,:),0); W(d,:) = W(d,:) + mu*(dbeta-S(d,:));
        M(d,:) = M(d,:) + mu*(dh - dX(d,:).*beta(d,:));
    end
    
    q = q + mu*sum(h,2); O = O + mu*(C-beta);
    Ln = max(abs(C),[],2); obj(iter) = 0.5*sum((y-sum(h,1)').^2) + lambda*sum(Ln);
    Ln1 = max(abs(beta),[],2);
    if mod(iter,100)==0
        disp(['    Slop:    Iter ' num2str(iter) ' Obj ' num2str(obj(iter))]);
        figure(9); subplot(2,1,1); plot(1:n,y,'ko',1:n,sum(h,1)','r+'); title('Response');
        subplot(2,1,2); plot(1:p,Ln,'r.',1:p,Ln1,'bo'); title('L\infty norm'); drawnow;
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

%%
clc; clear all; close all; format long; randn('state',0); rand('state',0);
p = 40; n = 200; k = 5; lambda = 10*sqrt(n); mu = 100;
X = randn(p,n);
Q =[1.0000    0.8000    0.6000    0.4000    0.2000
    0.8000    1.0000    0.7000    0.5000    0.3000
    0.6000    0.7000    1.0000    0.2000    0.1000
    0.4000    0.5000    0.2000    1.0000    0.1000
    0.2000    0.3000    0.1000    0.1000    1.0000]; 
% Q = eye(k);
ord = randperm(p); J = ones(p,1)==0; J(ord(1:k)) = 1==1;
y = sum(X(J,:).*(Q*X(J,:)),1)';
[beta,z,obj,Ln] = Piecewise_Slope(X,y-mean(y),lambda,mu,2000); Ln

