% SCAM_QP run
% X is Gaussian with covariance
% f is quadratic non-additive

p = 4;
n = 3000;
noise = 0.05;

% generate data
Z = randn(n,p);
tmp = randn(p,p);
covar = tmp'*tmp/(2*p) + eye(p)/2;

%covar = eye(p);
%covar = covar/(2*p);

X = Z * covar;

nout = sum(sum(abs(X) > 2));
X(abs(X) > 2) = rand(nout,1)-0.5;

% set up output
tmp_mat = randn(p,p);
fx_matrix = tmp_mat'*tmp_mat + eye(p);
fx_matrix(p,:) = 0;
fx_matrix(:,p) = 0;

fx_matrix = fx_matrix/(2*p);


fx = diag(X*fx_matrix*X');
y = fx + noise*randn(n,1);

meany = mean(y);
y = y - mean(y);

[beta,z,obj] = SCAM_QP(X', y, 0, 10, 1e-7);

for j = 1:p
    zgold = fx_matrix(j,j)*X(:,j).^2;
    
    [x_ord, ixs] = sort(X(:,j));
    
    figure;
    hold on;
    plot(x_ord, z(j,ixs), 'r');
    plot(x_ord, zgold(ixs), 'b');
    hold off;
    
end
