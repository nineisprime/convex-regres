%-------------------%
%
% Independent Script, uses only "lassopos"
%
%

Xs = randn(300,1);
[Xd, ord] = sort(Xs', 'ascend');
n = length(Xs);
xcol = reshape(Xs,n,1);

Xdrow = reshape(Xd(1:(n-1)), 1, n-1);

% form delta matrix
pieceI = repmat(xcol, 1, n-1);
pieceJ = repmat(Xdrow, n, 1);

Delta_matrix = max(pieceI - pieceJ, 0);
Delta2_matrix = Delta_matrix - (1/n)*repmat(sum(Delta_matrix, 1), n, 1);

X = Delta2_matrix;

betastar = [ones(1,4), zeros(1,299-4)];
betastar(20) = 1;
y = X*betastar';

beta = lassopos(X'*X, X'*y, 0, 300);
[beta,betastar']

