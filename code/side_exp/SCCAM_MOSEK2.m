function [beta1,beta2,z1,z2] = SCCAM_MOSEK2(y,x,lambda)
% Convex and Concave Regression
n = length(x); [Xd,ord] = sort(x); dX = Xd(2:n)-Xd(1:n-1);

Hi = [1:2*n 1:n n+1:2*n 2*n+1:4*n]; 
Hj = [1:2*n n+1:2*n 1:n 2*n+1:4*n];
Hs = [1/n*ones(1,2*n) -1/n*ones(1,n) -1/n*ones(1,n) 10^-6*ones(1,2*n)];
H = sparse(Hi,Hj,Hs,4*n,4*n); f = [-y/n; y/n; zeros(2*n-2,1); lambda; lambda];

Ai = []; Aj = []; As = []; row = 0; ndx = 0;
rg = (ndx+1):(ndx+2*(n-2)); ndx = ndx + 2*(n-2);
As(rg) = [ones(1,n-2) -ones(1,n-2)];
Ai(rg) = [1:(n-2) 1:(n-2)] + row; row = row + (n-2); 
Aj(rg) = [(2*n+1):(3*n-2) (2*n+2):(3*n-1)];

rg = (ndx+1):(ndx+4*(n-1)); ndx = ndx + 4*(n-1);
As(rg) = [ones(1,n-1) -ones(1,n-1) -ones(1,n-1) -ones(1,n-1)];
Ai(rg) = [1:(n-1) 1:(n-1) n:(2*n-2) n:(2*n-2)] + row; row = row + 2*(n-1); 
Aj(rg) = [(2*n+1):(3*n-1) (4*n-1)*ones(1,n-1) (2*n+1):(3*n-1) (4*n-1)*ones(1,n-1)];

rg = (ndx+1):(ndx+2*(n-2)); ndx = ndx + 2*(n-2);
As(rg) = [ones(1,n-2) -ones(1,n-2)];
Ai(rg) = [1:(n-2) 1:(n-2)] + row; row = row + (n-2); 
Aj(rg) = [3*n:(4*n-3) (3*n+1):(4*n-2)];

rg = (ndx+1):(ndx+4*(n-1)); ndx = ndx + 4*(n-1);
As(rg) = [ones(1,n-1) -ones(1,n-1) -ones(1,n-1) -ones(1,n-1)];
Ai(rg) = [1:(n-1) 1:(n-1) n:(2*n-2) n:(2*n-2)] + row; row = row + 2*(n-1); 
Aj(rg) = [3*n:(4*n-2) 4*n*ones(1,n-1) 3*n:(4*n-2) 4*n*ones(1,n-1)];

A = sparse(Ai,Aj,As,row,4*n); b = zeros(row,1);

Bi = []; Bj = []; Bs = []; row = 0; ndx = 0;
rg = (ndx+1):(ndx+n); ndx = ndx + n;
Bs(rg) = ones(1,n);
Bi(rg) = ones(1,n) + row; row = row + 1; 
Bj(rg) = 1:n;

rg = (ndx+1):(ndx+3*(n-1)); ndx = ndx + 3*(n-1);
Bs(rg) = [-ones(1,n-1) ones(1,n-1) -dX];
Bi(rg) = [1:n-1 1:n-1 1:n-1] + row; row = row + n-1; 
Bj(rg) = [ord(1:n-1) ord(2:n) (2*n+1):(3*n-1)];

rg = (ndx+1):(ndx+n); ndx = ndx + n;
Bs(rg) = ones(1,n);
Bi(rg) = ones(1,n) + row; row = row + 1; 
Bj(rg) = (n+1):2*n;

rg = (ndx+1):(ndx+3*(n-1)); ndx = ndx + 3*(n-1);
Bs(rg) = [-ones(1,n-1) ones(1,n-1) -dX];
Bi(rg) = [1:n-1 1:n-1 1:n-1] + row; row = row + n-1; 
Bj(rg) = [n+ord(1:n-1) n+ord(2:n) 3*n:(4*n-2)];

B = sparse(Bi,Bj,Bs,row,4*n); c = zeros(row,1);

h = quadprog(H,f,A,b,B,c); % MOSEK call of quadprog
z1 = h(1:n); z2 = h(n+1:2*n);
beta1 = h((2*n+1):(3*n-1)); beta2 = h(3*n:(4*n-2));
return