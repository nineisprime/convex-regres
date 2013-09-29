variable_selection_manual = 1;

p=200; k=10; s=5; n=200; sigma=0;
lambda = 0.001; mu=0.01; maxiter=300;

variable_selection_script;

figure;
plot(1:p, sum(abs(beta),2), '*r');
title('p=200, k=10, s=5, n=200, sigma=0');

%%%%%%

p=200; k=20; s=5; n=200; sigma=0;
lambda = 0.001; mu=0.01; maxiter=300;

variable_selection_script;

figure;
plot(1:p, sum(abs(beta),2), '*r');
title('p=200, k=20, s=5, n=200, sigma=0');

%%%%%%%%%%%%%

p=200; k=20; s=5; n=400; sigma=0;
lambda = 0.001; mu=0.01; maxiter=300;

variable_selection_script;

figure;
plot(1:p, sum(abs(beta),2), '*r');
title('p=200, k=20, s=5, n=400, sigma=0');

%%%%%%%%%%%

p=100; k=10; s=5; n=300; sigma=0.1;
lambda = 0.0001*n; mu=0.01; maxiter=500;

variable_selection_script;

figure;
plot(1:p, sum(abs(beta),2), '*r');
%title('p=100, k=10, s=5, n=300, sigma=0.1');

[sum(beta.^2,2),sum(beta_stars.^2,2)]

