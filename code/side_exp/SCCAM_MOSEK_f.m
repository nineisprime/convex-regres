function [z,z1,z2, Ln] = SCCAM_MOSEK_f(y,x,lambda)

n = length(x); [Xd,ord] = sort(x'); dX = Xd(2:n)-Xd(1:n-1);

yord = y(ord)';

%problem is that MOSEK fails constantly
%3*n + 1 variables

Hi = 1:3*n; 
Hj = 1:3*n; 
Hs = [ones(1,n)/n, 1e-6*ones(1,n), 1e-6*ones(1,n)];

H = sparse(Hi,Hj,Hs,3*n+1,3*n+1); f = [-yord/n; zeros(2*n, 1); lambda];

%put in Linf norm constraints
inf_mat1 = [eye(n), zeros(n), zeros(n), -ones(n,1)];
inf_mat2 = [-eye(n), zeros(n), zeros(n), -ones(n,1)];
A = sparse([inf_mat1; inf_mat2]);
% 
% inf_mat1 = [zeros(n), eye(n), zeros(n), -ones(n,1)];
% inf_mat2 = [zeros(n), -eye(n), zeros(n), -ones(n,1)];
 %inf_mat3 = [zeros(n), zeros(n), eye(n), -ones(n,1)];
 %inf_mat4 = [zeros(n), zeros(n), -eye(n), -ones(n,1)];
 %A = [inf_mat1; inf_mat2; inf_mat3; inf_mat4];
 %A = sparse(A);

%put in second-derivative constraints
adjmat = diag(1./dX,1) + diag(1./dX,-1);
degree_mat = diag(sum(adjmat,1));
laplace_mat = degree_mat - adjmat;
laplace_mat = laplace_mat(2:(n-1), :);

convex_mat1 = [zeros(n-2, n), laplace_mat, zeros(n-2,n), zeros(n-2,1)];
convex_mat1 = sparse(convex_mat1);
convex_mat2 = [zeros(n-2, n), zeros(n-2,n), laplace_mat, zeros(n-2,1)];
convex_mat2 = sparse(convex_mat2);

A = [A; convex_mat1; convex_mat2];
row = size(A,1);
b = zeros(row,1);

%put in zero-mean equality constraint
ave_mat1 = [zeros(1,n), ones(1,n), zeros(1,n), 0];
ave_mat2 = [zeros(1,n), zeros(1,n), ones(1,n), 0];
B = sparse([ave_mat1; ave_mat2]);

sum_mat = [eye(n), -eye(n), eye(n), zeros(n,1)];
sum_mat = sparse(sum_mat);

B = [B; sum_mat]; c = zeros(2+n,1); 

%put in summation equality

l = [-Inf*ones(3*n, 1); 0];

tstart=tic;
[h, fval, exitflag] = quadprog(H,f,A,b,B,c,l); % MOSEK call of quadprog

if (exitflag == -1)
    'MOSEK failed'
end

toc(tstart)
z = h(1:n); Ln = h(3*n+1);  
z1 = h((n+1):2*n);
z2 = h((2*n+1):3*n);


[~,ordrev] = sort(ord);
z = z(ordrev);
z1 = z1(ordrev);
z2 = z2(ordrev);




