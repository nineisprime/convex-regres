function [beta1,beta2,z1,z2] = SCCAM_MOSEK1_extend(y,x,lambda)
n = length(x); [Xd,ord] = sort(x); dX = Xd(2:n)-Xd(1:n-1);

%Hi = [1:n 1:n n+1:2*n 2*n+1:4*n-1]; 
%Hj = [1:n n+1:2*n 1:n 2*n+1:4*n-1];
Hi = 1:5*n-1;
Hj = 1:5*n-1;
Hs = [zeros(1,2*n) 10^-6*ones(1,2*n-1) ones(1,n)/n];
H = sparse(Hi,Hj,Hs,5*n-1,5*n-1); 

f = [zeros(2*n,1); zeros(2*n-2,1); lambda; -y/n];

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
Aj(rg) = [3*n:(4*n-2) (4*n-1)*ones(1,n-1) 3*n:(4*n-2) (4*n-1)*ones(1,n-1)];

A = sparse(Ai,Aj,As,row,5*n-1); b = zeros(row,1);

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

B = sparse(Bi,Bj,Bs,row,5*n-1); 

sum_mat = [-eye(n), eye(n), zeros(n, 2*n-1), eye(n)];
sum_mat = sparse(sum_mat);

B = [B; sum_mat];

row = size(B,1);
c = zeros(row,1);



%tstart = tic;
[h, fval, exitflag] = quadprog(H,f,A,b,B,c); % MOSEK call of quadprog
%toc(tstart)

if (exitflag == -1)
    'MOSEK FAILED'
end

z1 = h(1:n); z2 = h(n+1:2*n);
beta1 = h((2*n+1):(3*n-1)); beta2 = h(3*n:(4*n-2));
return
