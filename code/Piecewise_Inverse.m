function [ialpha,igamma,is] = Piecewise_Inverse(alpha,gamma,s)
% inv(diag(alpha)-s*gamma*gamma') = diag(ialpha)-is*igamma*igamma'.
ialpha = 1./alpha; igamma = gamma./alpha; is = - 1/(1/s-gamma'*igamma);
return