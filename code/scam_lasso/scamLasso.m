% OUTPUT:   beta is p--by--(n-1)
%           z is p--by--n
%
% NOTE: the following code is copied entirely from "SCAM_QP"
function [beta,z,obj,Ln] = scamLasso(X,y,lambda,maxit,tol)
[p,n] = size(X); z = zeros(p,n); 
Ln = zeros(p,1); beta = zeros(p,n-1);

res = y - sum(z,1)'; iter = 0; 
obj = []; change = 1;
counter = 0;
while iter < maxit && change > tol
    iter = iter + 1; 
    
    for d = 1:p
        if (mod(d,100)==0)
            d
            counter
        end
        res = res + z(d,:)';
        [beta2,z2,Ln2] = scamLassoOne(res,X(d,:),lambda/2);
        if (Ln2 < 1e-6)
            counter = counter+1;
            continue
        end
        [beta(d,:),z(d,:),Ln(d)] = SCAM_MOSEK(res,X(d,:),lambda);
        
        res = res - z(d,:)'; 
    end

    obj(iter) = 0.5/n*sum(res.^2) + lambda*sum(Ln); 
    if iter>1, change = abs(obj(iter-1)-obj(iter))/obj(iter-1); end
    disp(['   scamLasso:    Iteration ' num2str(iter) ' Objective ' num2str(obj(iter))]);
    figure(9); subplot(2,1,1); plot(1:n,y,'ko',1:n,sum(z,1)','r+'); title('Response');
    subplot(2,1,2); plot(1:p,Ln,'r.'); title('L\infty norm'); drawnow;
end
return



function [beta,z,Ln] = scamLassoOne(y,x,lambda)
n = length(x);
xcol = reshape(x, n, 1);

[Xd, ord] = sort(x', 'ascend');

Xdrow = reshape(Xd(1:(n-1)), 1, n-1);

% form delta matrix
pieceI = repmat(xcol, 1, n-1);
pieceJ = repmat(Xdrow, n, 1);

Delta_matrix = max(pieceI - pieceJ, 0);
Delta2_matrix = Delta_matrix - (1/n)*repmat(sum(Delta_matrix, 1), n, 1);

% form column for lasso
lasso_x = Delta2_matrix;
lasso_y = y;

xtx = lasso_x'*lasso_x;
xty = lasso_x'*lasso_y;

% do lasso
dvec = lassopos(xtx, xty, lambda, n);

% get beta
beta = cumsum(dvec);

% get z
z = Delta2_matrix * dvec;

% get Ln
Ln = max(abs(beta));







function [beta,z,Ln] = SCAM_MOSEK(y,x,lambda)
n = length(x); [Xd,ord] = sort(x'); dX = Xd(2:n)-Xd(1:n-1);

Hi = 1:n; Hj = 1:n; Hs = ones(1,n)/n; 
H = sparse(Hi,Hj,Hs,2*n,2*n); f = [-y/n; zeros(n-1,1); lambda];

Ai = []; Aj = []; As = []; row = 0; ndx = 0;
valmat = [ones(n-2,1) -ones(n-2,1)];
rowmat = [(1:(n-2))' (1:(n-2))'] + row;
colmat = [((n+1):(2*n-2))' ((n+2):(2*n-1))'];
rg = (ndx+1):(ndx+2*(n-2)); row = row + (n-2); ndx = ndx + 2*(n-2);
Ai(rg) = rowmat(:); Aj(rg) = colmat(:); As(rg) = valmat(:);

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

h = quadprog(H,f,A,b,B,c); % MOSEK call of quadprog
z = h(1:n); beta = h((n+1):(2*n-1)); Ln = h(2*n);  
return
