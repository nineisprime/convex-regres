% INPUT: X is (p--by--n)

% OUTPUT:   beta is p--by--(n-1)
%           z is p--by--n
%
function [beta,z,obj,Ln] = SCAM_QP(X,y,lambda,maxit,tol)
[p,n] = size(X); z = zeros(p,n); Ln = zeros(p,1); beta = zeros(p,n-1);

z2 = zeros(p,n);
Ln2 = zeros(p,1);

res = y - sum(z,1)'; iter = 0; obj = []; change = 1;
while iter < maxit && change > tol
    iter = iter + 1; 
    
    for d = 1:p
        res = res + z(d,:)';
        
        % MOSEK_f is much slower
        % both sometimes errs
        % 
        
        %disp('new')
        [z(d,:), Ln(d)] = SCAM_MOSEK_f(res, X(d,:), lambda);
        %disp('old')
        %[beta(d,:),z(d,:),Ln(d)] = SCAM_MOSEK(res,X(d,:),lambda);

        %[Xd, ord] = sort(X(d,:));
        %'test'
        %diff_magn = norm(z(d,:)-z2(d,:))/norm(z(d,:));
        %disp(['diff:' num2str( diff_magn )])
        
        %if (diff_magn > 1)
        %    'err'
        %end
        
        res = res - z(d,:)'; 
    end

    obj(iter) = 0.5/n*sum(res.^2) + lambda*sum(Ln); 
    if iter>1, change = abs(obj(iter-1)-obj(iter))/obj(iter-1); end
    disp(['   SCAM_QP:    Iteration ' num2str(iter) ' Objective ' num2str(obj(iter))]);
    figure(9); subplot(2,1,1); plot(1:n,y,'ko',1:n,sum(z,1)','r+'); title('Response');
    subplot(2,1,2); plot(1:p,Ln,'r.'); title('L\infty norm'); drawnow;
end
return

function [z,Ln] = SCAM_MOSEK_f(y,x,lambda)
n = length(x); [Xd,ord] = sort(x'); dX = Xd(2:n)-Xd(1:n-1);

yord = y(ord);

Hi = 1:n; Hj = 1:n; Hs = ones(1,n)/n; 
H = sparse(Hi,Hj,Hs,n+1,n+1); f = [-yord/n; lambda];

Ai = []; Aj = []; As = []; row = 0; 

%put in Linf norm constraints
valmat = [ones(n, 1) -ones(n, 1) -ones(n,1) -ones(n,1)];
rowmat = [1:n 1:n (n+1):(2*n) (n+1):(2*n)];
colmat = [1:n (n+1)*ones(1, n) 1:n (n+1)*ones(1,n)];
rg = 1:(4*n);  row = row + 2*n; 
Ai(rg) = rowmat'; Aj(rg) = colmat'; As(rg) = valmat(:);
A = sparse(Ai,Aj,As,row,n+1); 

%put in second-derivative constraints
adjmat = diag(1./dX,1) + diag(1./dX,-1);
degree_mat = diag(sum(adjmat,1));
laplace_mat = degree_mat - adjmat;
laplace_mat = laplace_mat(2:(n-1), :);
laplace_mat = [laplace_mat, zeros(n-2,1)];
laplace_mat = sparse(laplace_mat);

A = [A; laplace_mat];
row = row + (n-2);
b = zeros(row,1);

%put in zero-mean equality constraint
Bi = []; Bj = []; Bs = []; row = 0; 
valmat = ones(n,1); 
rowmat = ones(n,1) + row;
colmat = (1:n)';
rg = 1:n; row = row + 1; 
Bi(rg) = rowmat(:); Bj(rg) = colmat(:); Bs(rg) = valmat(:);
B = sparse(Bi,Bj,Bs,row,n+1); c = zeros(row,1); 

tstart=tic;
[h, fval, exitflag] = quadprog(H,f,A,b,B,c); % MOSEK call of quadprog

if (exitflag == -1)
    'MOSEK failed'
end

toc(tstart)
z = h(1:n); Ln = h(n+1);  

[~,ordrev] = sort(ord);
z = z(ordrev);

return

function [beta,z,Ln] = SCAM_MOSEK(y,x,lambda)
n = length(x); [Xd,ord] = sort(x'); dX = Xd(2:n)-Xd(1:n-1);

Hi = 1:n; Hj = 1:n; Hs = ones(1,n)/n; 
H = sparse(Hi,Hj,Hs,2*n,2*n); f = [-y/n; zeros(n-1,1); lambda];

%convexity
Ai = []; Aj = []; As = []; row = 0; ndx = 0;
valmat = [ones(n-2,1) -ones(n-2,1)];
rowmat = [(1:(n-2))' (1:(n-2))'] + row;
colmat = [((n+1):(2*n-2))' ((n+2):(2*n-1))'];
rg = (ndx+1):(ndx+2*(n-2)); row = row + (n-2); ndx = ndx + 2*(n-2);
Ai(rg) = rowmat(:); Aj(rg) = colmat(:); As(rg) = valmat(:);

%sup-norm
valmat = [ones(n-1,1) -ones(n-1,1) -ones(n-1,1) -ones(n-1,1)];
rowmat = [(1:(n-1))' (1:(n-1))' (n:(2*n-2))' (n:(2*n-2))'] + row;
colmat = [((n+1):(2*n-1))' 2*n*ones(n-1,1) ((n+1):(2*n-1))' 2*n*ones(n-1,1)];
rg = (ndx+1):(ndx+4*(n-1)); row = row + 2*(n-1); ndx = ndx + 4*(n-1);
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
tstart = tic;
[h, fval, exitflag] = quadprog(H,f,A,b,B,c); % MOSEK call of quadprog
if (exitflag == -1)
    'MOSEK failed'
end
toc(tstart)
z = h(1:n); beta = h((n+1):(2*n-1)); Ln = h(2*n);  
return