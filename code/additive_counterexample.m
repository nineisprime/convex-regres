
x1 = (0:1:3e1)/3e1;
x2 = (0:1:3e1)/3e1;

Y = zeros(length(x1),length(x2));
for ii = 1:length(x1)
    for jj=1:length(x2)
        Y(ii,jj) = sin(x1(ii)*2*pi)*sin(x2(jj)*2*pi);
    end
end
figure;
mesh(x1,x2,Y);
colormap(copper)


x1 = 2*(0:1:3e1)/3e1-1;
x2 = (0:1:3e1)/3e1;

Y = zeros(length(x1),length(x2));
for ii = 1:length(x1)
    for jj=1:length(x2)
        Y(ii,jj) = x1(ii)*x2(jj);
    end
end
figure;
mesh(x1,x2,Y);
colormap(copper);

n=5000;
Xs = [rand(n,1), 2*rand(n,1)-1];
Ys = Xs(:,1).*Xs(:,2);
[beta,z,obj] = SCAM_QP(Xs', Ys, 0.001, 100, 1e-5);


n=5000;
Xs = [rand(n,1), rand(n,1)];
Ys = sin(Xs(:,1)*(2*pi)).*sin(Xs(:,2)*(2*pi));
[beta,z,obj] = SCAM_QP(Xs', Ys, 0.001, 100, 1e-5);

