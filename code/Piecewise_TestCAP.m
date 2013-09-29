function err = Piecewise_TestCAP(type,n)
randn('state',0); rand('state',0);
ntrain = 5; err = zeros(ntrain,3);
switch type
    case 1
        p = 5; Xt = randn(p,n);
        yt = ((Xt(1,:) + 0.5*Xt(2,:) + Xt(3,:)).^2 - Xt(4,:) + 0.25*Xt(5,:).^2)';
    case 2
        p = 10; 
        prob = gamrnd(ones(p,1),1); prob = prob./sum(prob); 
        Xt = randn(p,n); yt = exp(Xt'*prob);
    otherwise, return
end

for run = 1:ntrain
    switch type
        case 1
            X = randn(p,n);
            y = ((X(1,:) + 0.5*X(2,:) + X(3,:)).^2 - X(4,:) + 0.25*X(5,:).^2)' + randn(n,1);
        case 2
            X = randn(p,n); y = exp(X'*prob) + 0.1*randn(n,1);
        otherwise, return
    end

    [alpha1,beta1] = Piecewise_CAP(X',y); % Lauren A. Hannah's code
    y1 = max([ones(n,1),Xt']*[alpha1; beta1],[],2);
    err(run,1) = sum((y1-yt).^2)/n;

    [beta2,h2,obj2] = Piecewise_Convex(X,y,0,0.01,200);
    % figure(1); plot(obj2,'r.'); drawnow;
    alpha2 = h2' - sum(X.*beta2,1);
    y2 = max([ones(n,1),Xt']*[alpha2; beta2],[],2);
    err(run,2) = sum((y2-yt).^2)/n;

    [beta3,h3,obj3] = Piecewise_QP(X,y);
    alpha3 = h3' - sum(X.*beta3,1);
    y3 = max([ones(n,1),Xt']*[alpha3; beta3],[],2);
    err(run,3) = sum((y3-yt).^2)/n;

    disp(['Finished type=' num2str(type) ' n=' num2str(n) ' run=' num2str(run)]);
end
return