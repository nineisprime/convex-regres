clc; clear all; close all;

p = 60;
K0 = 7;
s = 3;
K = 10; 

n_ls = [50,100,200,300,400,500];
noise=0.05;

num_trials = 10;
record_ls = zeros(length(n_ls),1);

for ii=1:length(n_ls)
    n = n_ls(ii);
    num_succ = 0;
    
    for it=1:num_trials    
        beta_stars = randn(p,K0);
        beta_stars((s+1):p,:)=0;
        alpha_stars = randn(K0,1)*0.2;
        Xs = randn(n,p); [Ys,c0] = max(Xs*beta_stars + ones(n,1)*alpha_stars', [], 2);
        [c0,ndx] = sort(c0); Xs = Xs(ndx,:); Ys = Ys(ndx) + noise*randn(n,1);

        c = randsample(K,n,true); 
        lambda = max(noise*sqrt(s*log(p)/n)*n,0.001*n); mu = 1; maxit = 1000;

        beta = zeros(p,K);
        alpha = zeros(1,K);
        for iter = 1:30
            [beta,alpha,obj,h,c,Ln] = Piecewise_KPlanesADMM(Xs',Ys,c,lambda,mu,maxit,beta,alpha);

            disp(['sparsity level: ' num2str(sum(max(abs(beta),[],2)<1e-5))]);

            disp(['Iteration ' num2str(iter) ' Objective ' num2str(obj(end))]);
            figure(1); subplot(3,2,1); plot(max(abs(beta_stars),[],2),'b*'); subplot(3,2,2); plot(Ln,'r.');
            subplot(3,2,3); plot(c0,'b.'); subplot(3,2,4); plot(c,'r.');
            subplot(3,2,[5,6]); plot(1:n,Ys,'ko',1:n,h,'r+'); drawnow;
        end
        close all;
        
        has_false_positive = sum(max(abs(beta((s+1):p,:)),[],2)) > 1e-2;
        has_false_negative = sum(max(abs(beta(1:s,:)),[],2)) < 1e-2;
        
        if (~has_false_positive && ~has_false_negative)
            num_succ = num_succ + 1;
        end
        
    end
    record_ls(ii) = num_succ/num_trials;
end
    