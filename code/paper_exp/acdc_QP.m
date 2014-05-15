% 

function [beta,z,obj,Ln] = acdc_QP(X,y,lambda,maxit,tol)
[p,n] = size(X); z = zeros(p,n); beta = zeros(p,n-1);

Ln = zeros(p,1); 
Ln2 = zeros(p,1);

res = y - sum(z,1)'; iter = 0; obj = []; change = 1;
obj = zeros(maxit,1);

counter = [0,0];
%SCAM
while iter < maxit && change > tol
    iter = iter + 1;
    
    for d = 1:p
        res = res + z(d,:)';
        
        
        % screening one
        xcol = reshape(X(d,:), n, 1);
        [Xd, ord] = sort(X(d,:), 'ascend');
        Xdrow = reshape(Xd(1:(n-1)), 1, n-1);
        pieceI = repmat(xcol, 1, n-1);
        pieceJ = repmat(Xdrow, n, 1);
        Delta_matrix = max(pieceI - pieceJ, 0);
        Delta2_matrix = Delta_matrix - (1/n)*repmat(sum(Delta_matrix, 1), n, 1);

        res_vec = (0.5/n)*res'*Delta_matrix;
        
        lambda_vec1 = lambda*(1/n)*ones(1,n)*Delta_matrix;
        
        uvec = abs(res'); uvec = uvec/sum(uvec);
        lambda_vec2 = lambda*uvec*Delta_matrix;
        
        screen_one_flag = sum(lambda_vec1 >= res_vec) == n-1 || ...
            sum(lambda_vec2 >= res_vec) == n-1 ;
        
        
   
        
        if (screen_one_flag)
            disp('good')
            counter(1) = counter(1) +1;
            %[beta(d,:),z(d,:),Ln(d)] = SCAM_MOSEK_f(res,X(d,:),lambda);
            %if (max(abs(z(d,:))) > 1e-6)
            %    disp('err!!!!!!!!')
            %end
            
            
        else
            disp('bad')
            counter(2) = counter(2) + 1;
            [beta(d,:),z(d,:),Ln(d)] = SCAM_MOSEK_f(res,X(d,:),lambda);
        end
        
        
        res = res - z(d,:)';
    end
    
    obj(iter) = 0.5/n*sum(res.^2) + lambda*sum(Ln);
    
    if iter>1
        change = abs(obj(iter-1)-obj(iter))/obj(iter-1); 
    end
    
    if mod(iter,10)==0
        disp(['   acdc_QP:    Iteration ' num2str(iter) ...
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
    
    xcol = reshape(X(d,:), n, 1);
    [Xd, ord] = sort(X(d,:), 'ascend');
    Xdrow = reshape(Xd(1:(n-1)), 1, n-1);
    pieceI = repmat(xcol, 1, n-1);
    pieceJ = repmat(Xdrow, n, 1);
    Delta_matrix = max(pieceI - pieceJ, 0);
    Delta2_matrix = Delta_matrix - (1/n)*repmat(sum(Delta_matrix, 1), n, 1);
    
    res_vec = -(.5/n)*res'*Delta_matrix;
    
    lambda_vec1 = lambda*(1/n)*ones(1,n)*Delta_matrix;
    
    uvec = abs(res'); uvec = uvec/sum(uvec);
    lambda_vec2 = lambda*uvec*Delta_matrix;
    
    screen_one_flag = sum(lambda_vec1 >= res_vec) == n-1 || ...
        sum(lambda_vec2 >= res_vec) == n-1 ;
    
    if (screen_one_flag)
        disp('good')
        counter(1) = counter(1) +1;
        %[~,~,Ln2(d)] = SCAM_MOSEK_f(res,X(d,:),lambda);
        %if (max(abs(z(d,:))) > 1e-6)
        %    disp('err')
        %end
        
        
    else
        disp('bad')
        counter(2) = counter(2) + 1;
        [~,~,Ln2(d)] = SCAM_MOSEK_f(-res,X(d,:),lambda);
    end
    
    
    res = res - z(d,:)';
end

Ln = max(Ln, Ln2);

disp(counter);

return

