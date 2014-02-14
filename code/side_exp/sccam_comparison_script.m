noise = 0.5;

n=200;

x = randn(n,1);
fx = sin(x*4*pi);
y = fx + noise*randn(n,1);

lambda = 0.01;

tstart = tic;
[beta1_a, beta2_a, z1_a, z2_a] = SCCAM_MOSEK1_extend(y, x', lambda);
toc(tstart)

tstart = tic;
[beta1_b, beta2_b, z1_b, z2_b] = SCCAM_MOSEK1(y, x', lambda);
toc(tstart)

diff_z1 = norm(z1_a - z1_b, 'fro')^2;
diff_z2 = norm(z2_a - z2_b, 'fro')^2;

norm_diff_z1 = diff_z1/(norm(z1_a,'fro')*norm(z1_b,'fro'))
norm_diff_z2 = diff_z2/(norm(z2_a,'fro')*norm(z2_b,'fro'))

max_diff_z1 = max(abs(z1_a - z1_b))
max_diff_z2 = max(abs(z2_a - z2_b))