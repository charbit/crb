%=============== transformlocsensor.m =========================
% to get isotropy for the slowness vector
% do t4isoflag = 1;
% The model is:
%      tau = xsensors_m * slowness + noise
% the covariance of the estimator of the D components of
% the slowness vector is proportional to identity.
% No simulation, only computation of the inv(HH^T)
%
clear
addpath toolbox
%====== flag to verify
t4isoflag     = 1;
%=======================
stationnumber = '37';
cdload        = sprintf('load ../../sensorlocation/I%s.mat',...
    stationnumber);
eval(cdload)
xsensors_m   = xsensors_m.coordinates;
xsensors_m   = xsensors_m(:,1:2);
[M,d]        = size(xsensors_m);
xsensors_m   = xsensors_m-ones(M,1)*mean(xsensors_m,1);
if t4isoflag
[xsensorsN_m,idkeep] = transform2isotrop(xsensors_m,1.2);
else
    xsensorsN_m = xsensors_m;
end
combi        = M*(M-1)/2;
distances    = zeros(combi,1);
orientations = zeros(combi,1);

cp=0;
for i1=1:M-1
    for i2=i1+1:M
        cp=cp+1;
        distances(cp)=norm(xsensorsN_m(i1,1:2)-xsensorsN_m(i2,1:2));
        orientations(cp)=atan((xsensorsN_m(i1,2)-...
            xsensorsN_m(i2,2))/(xsensorsN_m(i1,1)-xsensorsN_m(i2,1)));
    end
end
HHT       = inv(xsensorsN_m'*xsensorsN_m);
listalpha = 0:2:360;
Lalpha    = length(listalpha);
listv     = 340;
Lv        = length(listv);
res       = zeros(3,Lalpha,Lv);
for ia=1:Lalpha
    alpha = listalpha(ia)*pi/180;
    for iv=1:Lv
        v = listv(iv);
        J = ([-sin(alpha)/v -cos(alpha)/v/v;...
            cos(alpha)/v -sin(alpha)/v/v]);
        covmu = (J \ HHT) / J' ;
        res(:,ia,iv)=[covmu(1,1);covmu(2,2);covmu(1,2)];
    end
end

figure(2)
if t4isoflag
subplot(231)
plot(xsensors_m(idkeep(1:M-d),1)/1000,...
    xsensors_m(idkeep(1:M-d),2)/1000,'ob','markerfacec','r')
hold on
plot(xsensors_m(idkeep(M-d+1:M),1)/1000,...
    xsensors_m(idkeep(M-d+1:M),2)/1000,'ok','markerfacec','y')
hold on
plot(max(sqrt(sum(xsensors_m .^2,2)))*exp(2j*pi*(0:10:360)/360)/1000,':')
hold off
axis('square')
set(gca,'xlim',1.2*[-1, 1],'ylim',1.2*[-1, 1])
set(gca,'fontname','times','fontsize',12)
% xlabel('km')
title(sprintf('sensor locations\n%i sensor pairs',combi))

subplot(234)
plot(xsensorsN_m(1:M-d,1)/1000,xsensorsN_m(1:M-d,2)/1000,'ob','markerfacec','r')
hold on
plot(xsensorsN_m(M-d+1:M,1)/1000,xsensorsN_m(M-d+1:M,2)/1000,'ok','markerfacec','y')
plot(max(sqrt(sum(xsensorsN_m .^2,2)))*exp(2j*pi*(0:10:360)/360)/1000,':')
hold off
axis('square')
set(gca,'xlim',1.2*[-1, 1],'ylim',1.2*[-1, 1])
set(gca,'fontname','times','fontsize',12)
xlabel('km')
title(sprintf('moving the \n2 extreme sensors'))
else
    subplot(131)
    plot(xsensorsN_m(:,1)/1000,xsensorsN_m(:,2)/1000,'ob','markerfacec','r')
    hold on
    plot(max(sqrt(sum(xsensorsN_m .^2,2)))*exp(2j*pi*(0:10:360)/360)/1000,':')
    hold off
    axis('square')
    set(gca,'xlim',1.2*[-1, 1],'ylim',1.2*[-1, 1])
    set(gca,'fontname','times','fontsize',12)
    xlabel('km')
end
RR = [0.5, 1, 1.5 2];
subplot(133)
plot(0.1*(180/pi)*sqrt(squeeze(res(1,:,1))) .* exp(1j*listalpha*pi/180));
axis('square')
set(gca,'fontname','times','fontsize',12)
for ii=1:4
    hold on
    plot(RR(ii)*exp(2j*pi*(0:10:360)/360),':')
end
hold off
grid on
title(sprintf('typical accuracy\nas function of azimut'))

subplot(232)
plot(sort(distances)/1000,'.')
set(gca,'xlim',[1 combi])
title('inter-distances - km')
grid on

subplot(235)
plot(sort(orientations)*180/pi,'.')
set(gca,'ylim',90*[-1, 1])
set(gca,'xlim',[1 combi])
title('orientations - degree')
grid on

HorizontalSize = 15;
VerticalSize   = 11;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
% set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
%         set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');


printdirectory  = ' ../../slideslastPresentation/';
printfile = sprintf('%sstationexampleT.eps',printdirectory);
fileprintepscmd = sprintf('print -depsc -loose %s',printfile);
fileeps2pdfcmd  = sprintf('!epstopdf %s',printfile);
filermcmd       = sprintf('!rm %ss',printfile);

% eval(fileprintepscmd)
% eval(fileeps2pdfcmd)
% eval(filermcmd)