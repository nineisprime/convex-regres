function SCCAM2(noise)
% addpath('/Users/Minhua/mosek/6/toolbox/r2009b'); 
% setenv('MOSEKLM_LICENSE_FILE','/Users/Minhua/mosek/6/licenses/mosek.lic')
clc; close all; format long; randn('state',0); rand('state',0);
p = 64; n = 500; k = 5; 
X = randn(p,n); ord = randperm(p); h = zeros(k,n); hm = zeros(1,k); hstd = zeros(1,k);
sfun = {inline('abs(x)'),inline('abs(x).^2'),inline('x.^3'),inline('-abs(x).^4'),inline('exp(x)')};
for d = 1:k
    h(d,:) = feval(sfun{d},X(ord(d),:)); hm(d) = mean(h(d,:));
    hstd(d) = std(h(d,:)); h(d,:) = (h(d,:)-hm(d))/hstd(d);
end
y = sum(h,1)' + noise*randn(n,1);

lambda = 3*sqrt(log(n*p)/n); maxit = 20; tol = 10^-6; s = ones(p,2)==1;
[beta1,beta2,h1,h2,Ln,obj] = SCCAM_QP2(X,y,s,lambda,maxit,tol);
s = Ln>10^-8; active = find(sum(s,2)>0); num = length(active);
[beta1,beta2,h1,h2,Ln,obj] = SCCAM_QP2(X(active,:),y,s(active,:),0,maxit,tol);

xd = -2:0.001:2; hd(:,:,1) = SCCAM_Eval(repmat(xd,[num,1]),X(active,:),beta1,h1);
hd(:,:,2) = SCCAM_Eval(repmat(xd,[num,1]),X(active,:),beta2,h2);
for d = 1:num
    ndx = find(ord(1:k)==active(d)); yd = feval(sfun{ndx},xd);
    figure(999); subplot(2,3,d); plot(xd,hd(d,:,1),'b-',xd,-hd(d,:,2),'g-',xd,hd(d,:,1)-hd(d,:,2),'r-',...
        xd,(yd-hm(ndx))/hstd(ndx),'k-','LineWidth',2);
end
subplot(2,3,6); set(gca,'FontSize',14); plot(0,0,'b-',0,0,'g-',0,0,'r-',0,0,'k-'); axis off;
legend('Convex','Concave','Convex+Concave','Truth');
return