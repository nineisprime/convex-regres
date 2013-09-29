%Delta should be (p--by--n)

function [beta,s] = convexNoiseless(Delta, v, gamma)

    n = size(Delta,2);
    p = size(Delta,1);

    f = [zeros(1,p), 0*ones(1,n), ones(1,p)];
    %f = [zeros(1,p), 0*ones(1,n), 0.01*ones(1,3), ones(1,(p-3))];
    
    Aeq = zeros(n, n+2*p);
    Aeq(:,1:p) = Delta';
    Aeq(:,(p+1):(p+n)) = -eye(n);
    %Aeq(:,(n+p+1):(n+2*p)) = zeros(n,p);
    
    beq = v;
    
    A = zeros(n+2*p, n+p+p);
    for ii=1:p
        A(ii,ii) = 1;
        A(ii,ii+n+p) = -1;
        
        A(ii+n+p, ii) = -1;
        A(ii+n+p, ii+n+p) = -1;
    end
    
    for ii=1:n
        A(ii+p, ii+p) = 1;
    end
    
    b = zeros(n+2*p,1);
    
    [soln] = linprog(f,A,b,Aeq,beq);
    
    beta = soln(1:p);
    s = soln((p+1):(p+n));
    