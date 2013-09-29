x1 = (-100:1:100)/100;
x2=  [-1,-.6,-.1,.4,1]

y1 = -.5*(1-x1.^2).^(1/2);
y2 = -.5*(1-x2.^2).^(1/2);

figure;
hold on;
plot(x1,y1,'LineWidth',2);
plot(x2,y2,'k','LineWidth',2);
plot(x2,y2,'o','MarkerFaceColor','black','MarkerSize',7);

set(gca,'XTick',[]);
set(gca,'YTick',[]);
axis off;