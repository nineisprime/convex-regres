
% beta(2:p) must be positive
function [beta] = lassopos(xtx, xty, lambda, n)
%xtx = X'*X;
%xty = X'*Y;
%n = size(Y,1);
%eta = 0.9;
p = size(xtx, 1);

betax = zeros(p,1);

%betax_prev = betax;
betay = betax;
convergence_threshold = 1e-9;
betax_prev = betax;
t = 1;

opts.disp = 0;
stepsize = 1/eigs(xtx,1,'lm',opts);
%stepsize=1;

maxiter = 5e2;

for i = 1:maxiter
    
    %-------%
    % for FISTA
    %betay = betax + (i-2)/(i+1) * (betax - betax_prev);
    %betax_prev = betax;
    %--------%
    
    %stepsize = 1/i;
    grad = xtx*betay - xty;
    
    
    %%------------%
    %%for constant stepsize
    betax = project_l1(betay - stepsize*grad, n*stepsize*lambda);
    %betax(2:p) = max(betax(2:p), 0);
    %%----------%
    
    %        old_fobj = betay'*xtx*betay - 2*betay'*xty;
    %
    %         while true
    %             betax = project_l1(betay - stepsize*grad, n*stepsize*lambda);
    %             %betax = betay - stepsize*grad;
    %
    %             betax(2:p) = max(betax(2:p), 0);
    %
    %             new_fobj = betax'*xtx*betax - 2*betax'*xty;
    %
    %             if (new_fobj < old_fobj + (betax - betay)'*grad + ...
    %                     (1/(2*stepsize))*sum((betax - betay).^2))
    %                 break
    %             end
    %
    %             stepsize = stepsize*eta;
    %         end
    
    
    
    if (norm(betax - betax_prev, 2) < convergence_threshold)
        beta = betax;
        return
    end
    
    % for ISTA
    betax_prev = betax;
    betay = betax;
    
    %t_prev = t;
    %t = (1 + sqrt(1 + 4*t^2))/2;
    
end
beta = betax;

function [p] = project_l1(beta,lambda)
p = wthresh(beta,'s',lambda);

% X = randn(20, 50);
% y = X*[1, zeros(49)];
%
% beta = lassopos(