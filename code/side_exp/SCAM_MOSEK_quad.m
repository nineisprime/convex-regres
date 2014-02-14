
% Estimates y = f - c x^2 where c >= 0
% 
% L is an upper bound on c
%
% x is (1--by--n) row
%

function [beta,z,c,Ln] = SCAM_MOSEK_quad(y,x,lambda,L)
n = length(x); [Xd,ord] = sort(x'); dX = Xd(2:n)-Xd(1:n-1);

Hi = (2*n):(3*n-1); Hj = (2*n):(3*n-1); Hs = ones(1,n)/n; 
H = sparse(Hi,Hj,Hs,3*n+1,3*n+1); f = [ zeros(2*n-1,1); -y/n; 0; lambda];
%H(3*n, 3*n) = 1e-5/n;

%variable ordering
% [f..., beta...., h,...., c, gamma]

%convexity, constraint involve only beta
Ai = []; Aj = []; As = []; row = 0; ndx = 0;
valmat = [ones(n-2,1) -ones(n-2,1)];
rowmat = [(1:(n-2))' (1:(n-2))'] + row;
colmat = [((n+1):(2*n-2))' ((n+2):(2*n-1))'];
rg = (ndx+1):(ndx+2*(n-2)); row = row + (n-2); ndx = ndx + 2*(n-2);
Ai(rg) = rowmat(:); Aj(rg) = colmat(:); As(rg) = valmat(:);

A = sparse(Ai,Aj,As,row,3*n+1);

%sup-norm
block1 = [zeros(n-1, n), eye(n-1), zeros(n-1,n), -2*x(1:n-1)', -ones(n-1,1)];
block2 = [zeros(n-1, n), -eye(n-1), zeros(n-1,n), 2*x(1:n-1)', -ones(n-1,1)]; 
A = [A; sparse([block1; block2])];

row = size(A,1);
b = zeros(row,1);


Bi = []; Bj = []; Bs = []; row = 0; ndx = 0;
%mean 0, applies to h now
valmat = ones(n,1); 
rowmat = ones(n,1) + row;
colmat = (2*n:3*n-1)';
rg = (ndx+1):(ndx+n); row = row + 1; ndx = ndx + n;
Bi(rg) = rowmat(:); Bj(rg) = colmat(:); Bs(rg) = valmat(:);

valmat = [-ones(n-1,1) ones(n-1,1) -dX]; 
rowmat = [(1:n-1)' (1:n-1)' (1:n-1)'] + row;
colmat = [ord(1:n-1) ord(2:n) ((n+1):(2*n-1))'];
rg = (ndx+1):(ndx+3*(n-1)); row = row + n-1; ndx = ndx + 3*(n-1);
Bi(rg) = rowmat(:); Bj(rg) = colmat(:); Bs(rg) = valmat(:);

B = sparse(Bi,Bj,Bs,row,3*n+1); 

% substitution equalities
block = [ -eye(n), zeros(n,n-1), eye(n), x'.^2, zeros(n,1) ];
B = [B; block];

row = size(B,1);
c = zeros(row,1);

% enforce constraint on maximum value of quadratic-coeff c
l = -Inf*ones(3*n+1,1);
l(3*n) = 0;
u = Inf*ones(3*n+1,1);
u(3*n) = L;


tstart = tic;
[h, fval, exitflag] = quadprog(H,f,A,b,B,c,l,u); % MOSEK call of quadprog
if (exitflag == -1)
    'MOSEK failed'
end
toc(tstart)
z = h(2*n:3*n-1); beta = h((n+1):(2*n-1)); Ln = h(3*n+1);  
c = h(3*n);
return