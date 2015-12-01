% I37, I55, I27, I07

load('I37.mat')
% just load a file Ixx.mat
Msensors = length(xsensors_m.coordinates);
distorient = zeros(Msensors*(Msensors-1)/2,4);
grav_center = mean(xsensors_m.coordinates(:,1:2),1);
xsensors_m_centred = xsensors_m.coordinates(:,1:2) - ones(Msensors,1)*grav_center;

cp=0;
for i1=1:Msensors-1
    for i2=i1+1:Msensors
        cp=cp+1;
        Z = (xsensors_m.coordinates(i1,1) - xsensors_m.coordinates(i2,1)) + ...
            1j*(xsensors_m.coordinates(i1,2) - xsensors_m.coordinates(i2,2));
        distorient(cp,:) = [abs(Z) angle(Z)*180/pi i1 i2];
    end
end
C = cp;
[sortdistorient, indsortdistance] = sort(distorient,1);

    figure(1)
subplot(121)
plot(xsensors_m_centred(:,1)/1000,xsensors_m_centred(:,2)/1000,'o','markerfacec','r')
hold on
plot( max(sqrt((xsensors_m_centred .^2)*ones(2,1)))*exp(2j*pi*(0:100)/100)/1000,':')
set(gca,'xlim',1.100*[-1 1],'ylim',1.100*[-1 1])
hold off
axis('square')
set(gca, 'fontname','times','fontsize',12)
xlabel('km', 'fontname','times','fontsize',12)

subplot(222)
plot(sortdistorient(:,1)/1000,'.')
set(gca, 'fontname','times','fontsize',12)
% ylabel('km', 'fontname','times','fontsize',12)
set(gca,'xlim',[1 C])
grid on
title('inter-distances', 'fontname','times','fontsize',12)

subplot(224)
plot(sortdistorient(:,2),'.')
set(gca, 'fontname','times','fontsize',12)
set(gca,'xlim',[1 C],'ylim',180*[-1 1],'ytick',[-90 0 90])
grid on
% ylabel('degree', 'fontname','times','fontsize',10)
title('orientations', 'fontname','times','fontsize',12)

HorizontalSize = 18;
VerticalSize   = 10;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
        set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');


print -depsc ../myjob/2taskCRB/slides/stationexample.eps
!epstopdf ../myjob/2taskCRB/slides/stationexample.eps
!rm ../myjob/2taskCRB/slides/stationexample.eps


