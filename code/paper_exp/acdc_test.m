n = 1000;
alpha = 0.5;
%v = 0.3;
v = 0;

p = 512; 
k = 4; 

lambda = 0.5*sqrt(1/n)*log(n*p); 
maxit = 20; tol = 10^-6;

nrun = 1; 

J = ones(p,nrun)==0; 
Ln = zeros(p,nrun);


if alpha<10^-16
    Q = eye(k); 
else
    Q = 0.5*eye(k); 
    B0 = rand(k,k) < alpha;
    for i = 1:k
        for j = 1:i
            if B0(i,j) && i~=j
                Q(i,j) = 0.5; 
            end; 
        end; 
    end; Q = Q + Q';
end

for run = 1:nrun
    %ord = randperm(p); 
    ord = 1:p;
    J(ord(1:k),run) = 1==1; 
    
    %X = randn(p,n); 
    X = mvnrnd(zeros(1,p), toeplitz(v.^(0:p-1)), n)';
    
    y = sum(X(J(:,run),:).*(Q*X(J(:,run),:)),1)' + randn(n,1);
    % main function call
    [beta,h,obj,Ln(:,run)] = acdc_QP(X, y-mean(y), lambda, maxit, tol);
    %[beta,h,obj,Ln(:,run)] = SCAM_QP(X,y-mean(y),lambda,maxit,tol);
    %disp(['version=' num2str(version) ' run=' num2str(run)]);
end

prob = 0; epsil = 10^-6;

nrun = size(Ln,2); suc = 0;
for run = 1:nrun
    if max(Ln(~J(:,run),run)) < epsil && min(Ln(J(:,run),run)) > epsil
        suc = suc + 1;
    end
end
prob = suc/nrun;
