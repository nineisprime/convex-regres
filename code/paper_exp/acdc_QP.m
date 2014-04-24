% 

function [beta,z,obj,Ln] = acdc_QP(X,y,lambda,maxit,tol)
[p,n] = size(X); z = zeros(p,n); beta = zeros(p,n-1);

Ln = zeros(p,1); 
Ln2 = zeros(p,1);

res = y - sum(z,1)'; iter = 0; obj = []; change = 1;
obj = zeros(maxit,1);

%SCAM
while iter < maxit && change > tol
    iter = iter + 1;
    
    for d = 1:p
        res = res + z(d,:)';
        [beta(d,:),z(d,:),Ln(d)] = SCAM_MOSEK_f(res,X(d,:),lambda);
        res = res - z(d,:)';
    end
    
    obj(iter) = 0.5/n*sum(res.^2) + lambda*sum(Ln);
    
    if iter>1
        change = abs(obj(iter-1)-obj(iter))/obj(iter-1); 
    end
    
    if mod(iter,10)==0
        disp(['   SCAM_QP:    Iteration ' num2str(iter) ...
              ' Objective ' num2str(obj(iter))]);
        
        figure(9); subplot(2,1,1); 
        plot(1:n,y,'ko',1:n,sum(z,1)','r+'); 
        title('Response');
        
        subplot(2,1,2); plot(1:p,Ln,'r.'); 
        title('L\infty norm'); drawnow;
    end
end

%postprocessing
for d=1:p
    if Ln(d) > 1e-5
        continue
    end
    res = res + z(d,:)';
    [~, ~, Ln2(d)] = SCAM_MOSEK_f(-res, X(d,:), lambda);
    res = res - z(d,:)';
end

Ln = max(Ln, Ln2);


return

