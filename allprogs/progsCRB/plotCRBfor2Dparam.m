clear
alpha = 0.9;
N     = 1000;
X0    = [300;20];
W     = randn(N,2);
CRB   = 0.1*[0.2518, 0.2288;...
      0.2288, 1.0834];
sqrtCRB = sqrtm(CRB);
X     = W*sqrtCRB + ones(N,1)*X0';
figure(1)
plot(X(:,1),X(:,2),'o','color',0.7*ones(3,1),...
    'markerfacec',0.7*ones(3,1),'markers',4)
hold on
[area] = ellipse(X0, CRB, alpha,'b');
plot(X0(1),X0(2),'o','markers',12,'markerfacec','r')
hold off
xlabel('velocity - m/s')
ylabel('azimuth - degree')
set(gca,'xlim',X0(1)+[-1 1],'ylim',X0(2)+[-1 1])
set(gca','fontname','times','fontsize',14)

HorizontalSize = 8;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

print -depsc -loose ../../slideslastPresentation/crbexample.eps



