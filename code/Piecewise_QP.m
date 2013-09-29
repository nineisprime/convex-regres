function [beta,h,val] = Piecewise_QP(X,y)
[p,n] = size(X); % With l2 regularization on beta.
H = zeros(n*(p+1),n*(p+1)); H((n*p+(0:n-1))*(n*(p+1))+n*p+(1:n)) = 1;
H = sparse(H); 
f = [zeros(n*p,1);-y];
A = zeros(n*(n-1),n*(p+1));
row = 0;
for i = 1:n
    for j = 1:n
        if j~=i
            row = row + 1;
            A(row,((i-1)*p+1):(i*p)) = (X(:,j)-X(:,i))';
            A(row,n*p+j) = -1; A(row,n*p+i) = 1;
        end
    end
end
b = zeros(n*(n-1),1);
options = optimset('Algorithm','interior-point-convex','Display','off',...
                   'MaxIter',1000,'TolFun',10^-16,'TolX',10^-16);
[x,val] = quadprog(H,f,A,b,[],[],[],[],[],options); val = val + 0.5*y'*y;
beta = reshape(x(1:n*p),[p,n]); h = x(n*p+1:n*(p+1));
return