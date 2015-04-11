% OUT:
%   beta -- "p--by--(n-1)"
%   z -- "p--by--n"
%   used by "acdc_boston.m"

% "z" summed together estimates y

% IN:
%  "pattern" -- "p--by--1" vector of +/- 1


function [beta,z,obj] = acdc_refit(X,y,pattern,maxit,tol)

[p,n] = size(X); 
z = zeros(p,n); 
beta = zeros(p,n-1);



res = y; 
iter = 0; 
change = 1e5;
obj = zeros(maxit,1);


while iter < maxit && change > tol
    iter = iter + 1;
    
    for d = 1:p
        res = res + z(d,:)';
        
        if pattern(d) == 1
            %convex component
            [beta(d,:),z(d,:),~] = backfit_mosek_f(res,X(d,:),0);
        else
            %concave component
            [beta_cave,z_cave,~] = backfit_mosek_f(-res,X(d,:),0);
            z(d,:) = -z_cave;
            beta(d,:) = -beta_cave;
        end
            
        res = res - z(d,:)';
    end
    
    obj(iter) = 0.5/n*sum(res.^2);
    
    if iter>1
        change = abs(obj(iter-1)-obj(iter))/obj(iter-1); 
    end
    
    if mod(iter,2)==0
        disp(['   acdc_QP:    Iteration ' num2str(iter) ...
              ' Objective ' num2str(obj(iter))]);
        
        figure(9);  
        plot(1:n,y,'ko',1:n,sum(z,1)','r+'); 
        title('Response');
        
    end
end


return

