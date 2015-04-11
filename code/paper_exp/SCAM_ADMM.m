% ** DEPRECATED ** I believe
%
% Minhua's code, purpose unknown and undecipherable

function [C,h,obj,Ln] = SCAM_ADMM(X,y,lambda,mu,maxit)

[p,n] = size(X); 
dX = zeros(p,n-1); 
ord = zeros(p,n); 
Ln = zeros(p,maxit);

for d = 1:p, 
    [Xd,ord(d,:)] = sort(X(d,:)); 
    dX(d,:) = Xd(2:n)-Xd(1:n-1); 
end

beta = zeros(p,n-1); 
M = zeros(p,n-1); 
S = zeros(p,n-2); 
W = zeros(p,n-2);
C = zeros(p,n-1); 
O = zeros(p,n-1); 
q = zeros(p,1); 
iter = 0; 
obj = [];
G = gallery('tridiag',-ones(1,n-2),[2 3*ones(1,n-3) 2],-ones(1,n-2));
A = gallery('tridiag',-ones(1,n-1),[1 2*ones(1,n-2) 1],-ones(1,n-1));

B = inv(A+ones(n,n)); 
B1 = zeros(n,n); 
B2 = eye(n);

for d = 1:p, 
    B1(ord(d,:),:) = B; 
    B1(:,ord(d,:)) = B1; 
    B2 = B2+1/(mu*n)*B1; 
end; 

B2 = inv(B2);

rg = (0:n-2)*(n-1)+(1:n-1); 
ndx = (ord-1)*p + repmat((1:p)',[1 n]);
ndx1 = (ord(:,1:n-1)-1)*p + repmat((1:p)',[1 n-1]);
ndx2 = (ord(:,2:n)-1)*p + repmat((1:p)',[1 n-1]);

while iter < maxit
    iter = iter + 1;
    for d = 1:p, 
        C(d,:) = Piecewise_Inf(beta(d,:)-1/mu*O(d,:),lambda/mu); 
    end
    
    h = zeros(p,n); 
    eta = 1/(mu*n)*repmat(y',[p 1])-1/mu*repmat(q,[1 n]);
    resp = dX.*beta-1/mu*M; 
    eta(ndx) = eta(ndx) + [zeros(p,1) resp] - [resp zeros(p,1)];
    h(ndx) = eta(ndx)*B; 
    vec = B2*(sum(h,1)'); 
    h(ndx) = h(ndx) - 1/(mu*n)*vec(ord)*B;
    
    dh = h(ndx2) - h(ndx1);
    gd = C+1/mu*O + [zeros(p,1) S-1/mu*W]-[S-1/mu*W zeros(p,1)] + dX.*(dh+1/mu*M);
    for d = 1:p, 
        Gd = G; 
        Gd(rg) = Gd(rg) + dX(d,:).^2; 
        beta(d,:) = Gd\(gd(d,:)');
    end
    dbeta = beta(:,2:n-1) - beta(:,1:n-2);
    S = max(dbeta+1/mu*W,0); 
    W = W + mu*(dbeta-S); 
    M = M + mu*(dh - dX.*beta);
    
    q = q + mu*sum(h,2); 
    O = O + mu*(C-beta);
    Ln = max(abs(C),[],2); 
    obj(iter) = 0.5/n*sum((y-sum(h,1)').^2) + lambda*sum(Ln);
    Ln1 = max(abs(beta),[],2);
    if mod(iter,100)==0
        disp(['    Slop:    Iter ' num2str(iter) 
            ' Obj ' num2str(obj(iter))]);
        figure(9); 
        subplot(2,1,1); 
        plot(1:n,y,'ko',1:n,sum(h,1)','r+'); 
        title('Response');
        subplot(2,1,2); 
        plot(1:p,Ln,'r.',1:p,Ln1,'bo'); 
        title('L\infty norm'); 
        drawnow;
    end
end
return

function z = Piecewise_Inf(x,v)
n = length(x); 
[xa,ndx] = sort(abs(x),'descend');
b = [xa(2:n) 0] - (cumsum(xa)-v)./(1:n); 
z = zeros(1,n);

if b(n) < 0
    pos = find(b<0,1); 
    zinf = (sum(abs(x(ndx(1:pos))))-v)/pos;
    z(ndx(1:pos)) = sign(x(ndx(1:pos))).*zinf; 
    z(ndx(pos+1:n)) = x(ndx(pos+1:n));
end
return

%% Testing case:
clc; close all; format long; 
randn('state',0); rand('state',0);

n = 1000; p = 100;

k = 5; 
lambda = 3*sqrt(log(n*p)/n);
maxit = 10; tol = 10^-6; 

J = ones(p,1)==0; 
Q = eye(k); 

X = randn(p,n); 
ord = randperm(p); 
J(ord(1:k)) = 1==1;

y = sum(X(J,:).*(Q*X(J,:)),1)' + randn(n,1);

t = cputime; [beta1,h1,obj1,Ln1] = SCAM_QP(X,y-mean(y),lambda,maxit,tol); t1 = cputime - t;
t = cputime; [beta2,h2,obj2,Ln2] = SCAM_ADMM(X,y-mean(y),lambda,0.5,2000); t2 = cputime - t;

iter1 = length(obj1); iter2 = length(obj2); 
figure(1); subplot(2,1,1); 
plot(t1/iter1*(1:iter1),obj1,'r.-',t2/iter2*(1:iter2),obj2,'k.-');
xlabel('Elapsed time (seconds)'); 
ylabel('Objective value'); 
legend('QP','ADMM');

subplot(2,1,2); 
plot(1:p,Ln1,'r.',1:p,Ln2,'k.'); 
xlabel('Feature dimension'); 
ylabel('|\beta|_\infty'); 
legend('QP','ADMM');
