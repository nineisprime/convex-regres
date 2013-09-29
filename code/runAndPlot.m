function [out] = runAndPlot(params)

n=params.n;   p=params.p;   noise=params.noise;
s=params.s;
gamma = params.gamma;
fig_name = params.fig_name;

use_quad=params.use_quad;   use_lin=params.use_lin;
use_cube=params.use_cube;   use_logexp=params.use_logexp;

maxit = 3000;
mu = .1;

lambda = noise*sqrt(s*log(p)/n)*n; 
lambda = max(lambda,1e-3*n);

X = randn(p,n)/2;

K = 4;
beta_star = randn(p,K);
beta_star((s+1):p,:)=0;
y_lin = max(X'*beta_star,[],2);

tmp = randn(s,s);
A = tmp*tmp' + 0.5*eye(s);
A = A/norm(A,2);
y_quad = diag(X(1:s,:)'*A*X(1:s,:));

y_cube = sum(abs(X(1:s,:)).^3,1)';

beta_exp = rand(s,1)+0.2;
y_logexp = log(sum(exp((beta_exp*ones(1,n)).*X(1:s,:)),1))';

y = use_lin*y_lin + use_quad*y_quad + ...
    use_cube*y_cube + use_logexp*y_logexp + noise*randn(n,1);

[beta1,D,h1,obj1] = Quadratic_Infinity(X,y,lambda,gamma,mu,maxit);
[beta2,h2,obj2] = Piecewise_Infinity(X,y,lambda,mu,maxit);

fig_handle = figure;
% plot things here!
C = [beta1,D];
weights1 = max(abs(C),[],2);
weights2 = max(abs(beta2),[],2);

weights1 = min(weights1,0.5);
weights2 = min(weights2,0.5);

hold on;
plot(1:p, weights1,'or','MarkerSize',7);
plot(1:p, weights2, '*b', 'MarkerSize',7);

mean1 = [mean(weights1(1:s))*ones(s,1); mean(weights1(s+1:p))*ones(p-s,1)];
mean2 = [mean(weights2(1:s))*ones(s,1); mean(weights2(s+1:p))*ones(p-s,1)];

plot(1:p, mean1, '--r', 'LineWidth',1);
plot(1:p, mean2, '--b', 'LineWidth',1);

title(strcat('n=',num2str(n),' p=',num2str(p),' s=',num2str(s),' noise=',num2str(noise),...
    ' gamma=',num2str(gamma),' useQuad=',num2str(use_quad),...
    ' useLin=',num2str(use_lin),...
    ' useCube=',num2str(use_cube),...
    ' useLogexp=',num2str(use_logexp)));
xlabel('coordinate indices');
ylabel('coordinate weights');
legend('Quadratic','Piecewise Lin');
set(gca,'FontSize',14);
print(fig_handle, '-depsc', strcat('figs\Figure',fig_name));

out = 1;