n = 40;
p = 200;

s = 1;

X = randn(p,n);
%X = X./(sqrt(sum(X.^2,2))*ones(1,n));  

A = eye(s);
%y = (1/2)*diag(X(1:s,:)'*A*X(1:s,:));
%y = randn(n,1);
%beta = zeros(p,1);
%beta(1:s) = 0.5;
%noise = rand(n,1)*2;
%noise(1)=0;
%y = X'*beta + noise;
y = X(1,:)';

Delta = zeros(p,n);

%v = zeros(n,1);

%for ii=1:n
%    Delta(:,ii) = X(:,ii) - X(:,1);
%    v(ii) = y(ii) - y(1);
%end

%[beta,S] = convexNoiseless(Delta,v,0);

[beta,S] = convexNoiseless(X,y,0);

norm(beta,1)
norm(X(1:s,1),1)

plot(abs(beta),'r.')

