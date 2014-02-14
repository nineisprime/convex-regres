%function [z,Ln,obj] = SCCAM_QP1(X,y,lambda,maxit,tol)
function [beta1,beta2,z1,z2,Ln,obj] = SCCAM_QP1_beta(X,y,lambda,maxit,tol)
[p,n] = size(X); 

z1 = zeros(p,n); z2 = zeros(p,n); 
z = zeros(p,n);

beta1 = zeros(p,n-1); beta2 = zeros(p,n-1); 
Ln = zeros(p,1); 

res = y - sum(z1-z2,1)'; 
%res = y' - sum(z,1);

iter = 0; obj = []; change = 1;
while iter < maxit && change > tol
    iter = iter + 1;
    
    for d = 1:p
        res = res + (z1(d,:)-z2(d,:))';
        %res = res + z(d,:);
        tstart = tic;
        [beta1(d,:),beta2(d,:),z1(d,:),z2(d,:)] = SCCAM_MOSEK1_extend(res,X(d,:),lambda);
        %[z(d,:), z1, z2, ~] = SCCAM_MOSEK_f(res, X(d,:), lambda);
        toc(tstart)
        
        %res = res - z(d,:);
        res = res - (z1(d,:)-z2(d,:))';
    end
    
    Ln = max(abs([beta1 beta2]),[],2); 
    %Ln = max(abs(z), [], 2);
    obj(iter) = 0.5/n*sum(res.^2) + lambda*sum(Ln);
    if iter>1, change = abs(obj(iter-1)-obj(iter))/obj(iter-1); end
    if mod(iter,1)==0
        disp(['   SCCAM_QP1:    Iteration ' num2str(iter) ' Objective ' num2str(obj(iter))]);
        figure(9); subplot(2,1,1); 
        plot(1:n,y,'ko',1:n,sum(z1-z2,1)','r+'); 
        %plot(1:n,y,'ko',1:n,sum(z,1)','r+'); 
        title('Response');
        subplot(2,1,2); plot(1:p,Ln,'r.','markersize',10); title('L\infty norm'); drawnow;
    end
end
return
