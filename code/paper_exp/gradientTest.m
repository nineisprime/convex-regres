
n=41;
X = rand(n,1);

xcol = X;
[Xd, ord] = sort(X, 'ascend');
Xdrow = reshape(Xd(1:(n-1)), 1, n-1);
pieceI = repmat(xcol, 1, n-1);
pieceJ = repmat(Xdrow, n, 1);
Delta_matrix = max(pieceI - pieceJ, 0);
Delta_matrix = Delta_matrix(ord,:);
Delta2_matrix = Delta_matrix - (1/n)*repmat(sum(Delta_matrix, 1), n, 1);


dtd = Delta2_matrix'*Delta2_matrix;

rhat = randn(n,1);
rhat = rhat-mean(rhat);
tmp = (1/n)*Delta_matrix'*rhat;



%u = u*abs(tmp(1))/tmp2(1);
%tmp2 = (5/sqrt(n))*Delta_matrix'*u;
%if (tmp2(1) > tmp(1))
%    u(n) = -1*u(n-1);
%end
%tmp2 = (5/sqrt(n))*Delta_matrix'*u;

% begin mosek program
f = [zeros(n,1); ones(n,1)];

nvars = 2*n;
A = zeros(n-2 + 2*n, nvars);
A(1:(n-2), 1:n) = -(7/sqrt(n))*Delta_matrix(1:n, 2:(n-1))';
b = -[(1/n)*Delta_matrix(1:n, 2:(n-1))'*rhat; zeros(2*n,1)];

A((n-2)+(1:n), 1:nvars) = [eye(n), -eye(n)];
A((n-2+n)+(1:n), 1:nvars) = [-eye(n),-eye(n)];

B = zeros(1, nvars);
B(1, 1:n) = (7/sqrt(n))*Delta_matrix(1:n, 1)';
c = (1/n)*Delta_matrix(1:n, 1)'*rhat;

l = -Inf*ones(nvars,1);
u = [zeros(n-1,1); Inf*ones(n+1,1)];
l(1:floor(n/2)) = 0;
l((floor(n/2)+2):n)=0;


[x, fval, exitflag] = linprog(f,A,b,B,c,l,u);

u = x(1:n);
l1tot = sum(abs(u));

tmp2 = (7/sqrt(n))*Delta_matrix'*u;
[tmp, tmp2]
l1tot
plot(u)


M = tril(ones(n-1, n));
M = M(1:(n-1), fliplr(1:n));
%M*rhat

