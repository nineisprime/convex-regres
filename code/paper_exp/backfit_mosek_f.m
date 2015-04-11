% IN: 
%       y -- "n" length
%       x -- "n" length
%       lambda -- scalar
%
% OUT:
%       beta -- "n-1" length derivatives
%       z -- "n" length fitted values
%       Ln -- scalar for function Linf value
%
% performs BACKFIT with L1 penalty on function value
% 
% minimize over {f, beta, gamma}
%


function [beta,z,Ln] = backfit_mosek_f(y,x,lambda)
n = length(x); 
[Xd,ord] = sort(x'); 
dX = Xd(2:n)-Xd(1:n-1);

Hi = 1:2*n; 
Hj = 1:2*n; 
Hs = [ones(1,n)/n,1e-10*ones(1,n)];

% Use 1e-10 instead of 0 to make this QP strictly convex.
H = sparse(Hi,Hj,Hs,2*n,2*n); 
f = [-y/n; zeros(n-1,1); lambda];
Ai = []; 
Aj = []; 
As = []; 
row = 0; 
nnz_in_A = 0;

%variable ordering:
% function values: 1:n
% beta: (n+1):(2*n-1)
% gamma: 2*n
% A is (num constraint) * (num variables)

%% inequality constraints
% beta increasing constraints
% (1, n+1) -> 1    (1, n+2) -> -1
% (2, n+2) -> 1    (2, n+3) -> -1
% ...
%
% A is  1, -1
%           1, -1 ...
new_nnz = 2*(n-2);
new_constraints = n-2;

valmat = [ones(n-2,1) -ones(n-2,1)];
rowmat = [(1:(n-2))' (1:(n-2))'] + row;
colmat = [((n+1):(2*n-2))' ((n+2):(2*n-1))'];


rg = (nnz_in_A + 1):(nnz_in_A + new_nnz); 

row = row + new_constraints;
nnz_in_A = nnz_in_A + new_nnz;
Ai(rg) = rowmat(:); 
Aj(rg) = colmat(:); 
As(rg) = valmat(:);



% Linf penalty constraints
% first n variables must be <= gamma (2*n-th variable)
%                        and >= -gamma
    
% (1, 1) -> 1     (1, 2*n) -> -1
% (2, 2) -> 1     (2, 2*n) -> -1
% ... 
% (1, 1) -> -1    (1, 2*n) -> -1
% (2, 2) -> -1    (2, 2*n) -> -1

new_nnz = 4*n;   
new_constraints = 2*n;

valmat = [ones(n,1)  -ones(n,1) ...
          -ones(n,1)  -ones(n,1)];
rowmat = [(1:n)' (1:n)' ((n+1):(2*n))' ...
          ((n+1):(2*n))'] + row;
colmat = [(1:n)' 2*n*ones(n, 1) ...
          (1:n)' 2*n*ones(n, 1)];
   
rg = (nnz_in_A + 1):(nnz_in_A + new_nnz); 
row = row + new_constraints;
Ai(rg) = rowmat(:); 
Aj(rg) = colmat(:); 
As(rg) = valmat(:);

A = sparse(Ai,Aj,As,row,2*n); b = zeros(row,1);


%% Equality constraints
Bi = []; 
Bj = []; 
Bs = []; 
row = 0; 
nnz_in_B = 0;


% centering constraint
% (1, 1) -> 1, (1, 2) -> 1, ...

new_nnz = n;
new_constraints = 1;

valmat = ones(n,1);
rowmat = ones(n,1) + row;
colmat = (1:n)';

rg = (nnz_in_B + 1):(nnz_in_B + new_nnz); 
row = row + new_constraints; 
nnz_in_B = nnz_in_B + new_nnz;
Bi(rg) = rowmat(:); 
Bj(rg) = colmat(:); 
Bs(rg) = valmat(:);


% linear interpolation constraints
% (1, ord(1)) -> -1   (1, ord(2)) -> 1   (1, n+1) -> -dX(1)
% (2, ord(2)) -> -1   (2, ord(3)) -> 1   (2, n+2) -> -dx(2)

new_nnz = 3*(n-1);
new_constraints = n-1;

rowmat = [(1:n-1)' (1:n-1)' (1:n-1)'] + row;
colmat = [ord(1:n-1) ord(2:n) ((n+1):(2*n-1))'];
valmat = [-ones(n-1,1) ones(n-1,1) -dX];

rg = (nnz_in_B + 1):(nnz_in_B + new_nnz); 
row = row + new_constraints;

Bi(rg) = rowmat(:); 
Bj(rg) = colmat(:); 
Bs(rg) = valmat(:);

B = sparse(Bi,Bj,Bs,row,2*n); c = zeros(row,1);

h = quadprog(H,f,A,b,B,c); % MOSEK call of quadprog
z = h(1:n); beta = h((n+1):(2*n-1)); Ln = h(2*n);
return
