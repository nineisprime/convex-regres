
%
% performs SCAM_MOSEK with L1 penalty on function value
% 
% minimize over {f, beta, gamma}

function [beta,z,Ln] = SCAM_MOSEK_f(y,x,lambda)
n = length(x); [Xd,ord] = sort(x'); dX = Xd(2:n)-Xd(1:n-1);
Hi = 1:2*n; Hj = 1:2*n; Hs = [ones(1,n)/n,10^-6*ones(1,n)];

% Use 10^-6 instead of 0 to make this QP strictly convex.
H = sparse(Hi,Hj,Hs,2*n,2*n); f = [-y/n; zeros(n-1,1); lambda];
Ai = []; Aj = []; As = []; row = 0; ndx = 0;

%variable ordering:
% function values: 1:n
% beta: (n+1):(2*n-1)
% gamma: 2*n

valmat = [ones(n-2,1) -ones(n-2,1)];
rowmat = [(1:(n-2))' (1:(n-2))'] + row;
colmat = [((n+1):(2*n-2))' ((n+2):(2*n-1))'];
rg = (ndx+1):(ndx+2*(n-2)); row = row + (n-2); ndx = ndx + 2*(n-2);
Ai(rg) = rowmat(:); Aj(rg) = colmat(:); As(rg) = valmat(:);




% gamma inequalities
%valmat = [ones(n-1,1) -ones(n-1,1) -ones(n-1,1) -ones(n-1,1)];
valmat = [ones(n,1), -ones(n,1), -ones(n,1), -ones(n,1)];
%rowmat = [(1:(n-1))' (1:(n-1))' (n:(2*n-2))' (n:(2*n-2))'] + row;
rowmat = [(1:n)' (1:n)' ((n+1):(2*n))' ((n+1):(2*n))'] + row;
%colmat = [((n+1):(2*n-1))' 2*n*ones(n-1,1) ...
%          ((n+1):(2*n-1))' 2*n*ones(n-1,1)];
colmat = [(1:n)' 2*n*ones(n, 1) ...
          (1:n)' 2*n*ones(n, 1)];
%rg = (ndx+1):(ndx+4*(n-1)); row = row + 2*(n-1); ...
%      ndx = ndx + 4*(n-1);
rg = (ndx+1):(ndx+4*n); 
row = row + 2*n;
Ai(rg) = rowmat(:); Aj(rg) = colmat(:); As(rg) = valmat(:);

A = sparse(Ai,Aj,As,row,2*n); b = zeros(row,1);





Bi = []; Bj = []; Bs = []; row = 0; ndx = 0;
valmat = ones(n,1);
rowmat = ones(n,1) + row;
colmat = (1:n)';
rg = (ndx+1):(ndx+n); row = row + 1; ndx = ndx + n;
Bi(rg) = rowmat(:); Bj(rg) = colmat(:); Bs(rg) = valmat(:);

valmat = [-ones(n-1,1) ones(n-1,1) -dX];
rowmat = [(1:n-1)' (1:n-1)' (1:n-1)'] + row;
colmat = [ord(1:n-1) ord(2:n) ((n+1):(2*n-1))'];
rg = (ndx+1):(ndx+3*(n-1)); row = row + n-1; ndx = ndx + 3*(n-1);
Bi(rg) = rowmat(:); Bj(rg) = colmat(:); Bs(rg) = valmat(:);

B = sparse(Bi,Bj,Bs,row,2*n); c = zeros(row,1);

h = quadprog(H,f,A,b,B,c); % MOSEK call of quadprog
z = h(1:n); beta = h((n+1):(2*n-1)); Ln = h(2*n);
return
