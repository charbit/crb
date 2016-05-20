%======================= CRBwithLOC.m
clear

addpath toolbox

allcolors = ['g.';'m.';'r.';'k.';'b.';'rx';...
    'yx';'mx';'rx';'kx';'c.';'k.';'r.';'c.';'m.';...
    'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

listaz = 0; %linspace(0,0,1)';
Laz = length(listaz);
stdaz = zeros(Laz,1);
stdel = zeros(Laz,1);
stdvel = zeros(Laz,1);

aec.e_deg      = 20;
aec.c_mps      = 340;
Fs_Hz          = 20;
SNR_dB         = 0;
T_cal          = 100;
N              = fix(T_cal*Fs_Hz);
K              = fix(N/2);
frequency_Hz   = (1:K)'*Fs_Hz/N;
sigma2noise    = 10^(-SNR_dB/10);
modifloc       = 0;
alpha_coh      = 0;%5e-5;%9e-5;%0;%1.42e-4;%0.008;
magnEV         = 1.2;
Llistfactor    = 12;
listfactor     = linspace(500,2000,Llistfactor);
numfig         = 1;
numfig2        = 3;
printflag      = 0;

choice = 3;
switch choice
    case 5
        M               = 8;
        xsensors_m        = zeros(M,3);
        xsensors_m(:,1)   = (-M+1:2:M-1)'/M;
    case 1
        M               = 3;
        aux             = exp(2j*pi*(0:M-1)'/M);
        xsensor0(:,1)   = real(aux);
        xsensor0(:,2)   = imag(aux);
        xsensor0(:,3)   = 1*ones(M,1);
    case 2
        M               = 8;
        Mon2            = M/2;
        aux1            = exp(2j*pi*(0:Mon2-1)'/Mon2);
        aux2            = exp(2j*pi*(0:Mon2-1)'/Mon2+1j*pi/4)/2;
        xsensor0(:,1)   = real([aux1;aux2]);
        xsensor0(:,2)   = imag([aux1;aux2]);
        xsensor0(:,3)   = 0*ones(M,1);
    case 6
        stationnumber            = 37;
        comm = sprintf('load ../../sensorlocation/I%i.mat',stationnumber);
        eval(comm)
        xsensors_m               = xsensors_m.coordinates(:,1:2);
        [M,D]                    = size(xsensors_m);
        xsensors_centered_m      = xsensors_m-ones(M,1)*mean(xsensors_m,1);
        R0                       = sqrt(max(sum(xsensors_centered_m .^2,2)));
        xsensors_centered_norm_m = xsensors_centered_m/R0;
    case 3
        xsensors_m =  1000*[
            [-0.05997213,  0.194591122,   0.3911];
            [0.229169719,   0.083396195,   0.3921];
            [0.122158887,  -0.206822564,   0.3918];
            [-0.12375342,  -0.087843992,   0.3902];
            [-0.026664123,   0.015567290,   0.391];
            [0.919425013,   0.719431175,   0.3924];
            [0.183105453,  -1.103053672,   0.3831];
            [-1.2434694,   0.384734446,   0.3976]];
        [M,D] = size(xsensors_m);
        %
end

if modifloc
    xsensors_new_m = transform2isotrop(xsensors_centered_norm_m,magnEV);
else
    xsensors_new_m = xsensors_m ;%xsensors_centered_norm_m;
end
listlegs = cell(Llistfactor,1);
meanstdaz = zeros(Llistfactor,1);
meanstdvel = zeros(Llistfactor,1);
for ifactor=1:Llistfactor
    
    factor    = listfactor(ifactor);
    listlegs{ifactor} = sprintf('%s:%i m',allcolors(ifactor), factor);
    xsensors_fact_m   = factor * xsensors_new_m;
    sensordistance = zeros(M,M);
    for im1=1:M
        for im2=1:M
            sensordistance(im1,im2)=norm(xsensors_fact_m(im1,:)-...
                xsensors_fact_m(im2,:));
        end
    end
    C     = ones(K,M,M);
    for ik=1:K
        for im1=1:M-1
            for im2=im1+1:M
                C(ik,im1,im2)=exp(-alpha_coh*(sensordistance(im1,im2)^2*frequency_Hz(ik)^2 ));
                C(ik,im2,im1)=C(ik,im1,im2);
            end
        end
    end
    for iaz=1:Laz
        aec.a_deg = listaz(iaz);
        switch D
            case 3
                polparameters = aec;
            case 2
                polparameters.a_deg = aec.a_deg;
                polparameters.v_mps = aec.c_mps;
        end
        CRB          = evalCRBwithLOC(xsensors_m, sigma2noise, C, polparameters, frequency_Hz);
        stdaz(iaz)   = sqrt(CRB.av(1,1))*180/pi;
        stdvel(iaz)  = sqrt(CRB.av(2,2));
        
    end
    azasaz = stdaz .* exp(1j*pi*listaz/180);
    meanstdaz(ifactor) = mean(stdaz);
    velasaz = stdvel .* exp(1j*pi*listaz/180);
    
    meanstdvel(ifactor) = mean(stdvel);
    
    
    subplot(121)
    if not(ifactor==1)
        hold on
    end
    plot(azasaz,'.-','color',allcolors(ifactor))
    hold off
    
    subplot(122)
    if not(ifactor==1)
        hold on
    end
    plot(velasaz,'.-','color',allcolors(ifactor))
    hold off
    
end
%%
figure(numfig)

subplot(121)
set(gca,'fontname','times','fontsize',12)
% Mmax = 2;
% set(gca,'xlim',Mmax*[-1,1])
% set(gca,'ylim',Mmax*[-1,1])
% set(gca,'xtick',[0],'ytick',[0])
xlabel('azimuth');
axis('square')
grid on

%==
subplot(122)
set(gca,'fontname','times','fontsize',12)


% Mmax = 10;
% set(gca,'xlim',Mmax*[-1,1])
% set(gca,'ylim',Mmax*[-1,1])
% set(gca,'xtick',[0],'ytick',[0])
xlabel('hor. velocity');
axis('square')
grid on

% subplot(133)
% plot(xsensors_fact_m(:,1),xsensors_fact_m(:,2),'o')
% subplot(133)
% % Mmax = 1;
% % set(gca,'xlim',Mmax*[-1,1])
% % set(gca,'ylim',Mmax*[-1,1])
% axis('square')
% % set(gca,'xtick',[0],'ytick',[0])
% grid on
% % CRB.slowness
% % CRB.aec
% % sqrt(CRB.aec(1,1))*180/pi
% % sqrt(CRB.aec(2,2))*180/pi
% % sqrt(CRB.aec(3,3))
% %
figure(numfig)
HorizontalSize = 24;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);

set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
set(gcf, 'InvertHardCopy', 'off');



%%
figure(numfig2)
subplot(131)
plot(xsensors_new_m(:,1),xsensors_new_m(:,2),'o',...
    'markersize',12,'markerfacec','r')
axis('square')

hold on
plot(max(sqrt(sum(xsensors_new_m .^2,2)))*exp(2j*pi*(0:10:360)/360),':')
hold off
set(gca,'fontname','times','fontsize',12)
xlabel('m')
ylabel('m')

subplot(132)
plot(listfactor,meanstdaz,'.-')
set(gca,'fontname','times','fontsize',12)
xlabel('multiplying factor')
ylabel('azimuth STD mean')
grid on
title(sprintf('LOC-coeff = %4.2e',alpha_coh))

subplot(133)
plot(listfactor,meanstdvel,'.-')
set(gca,'fontname','times','fontsize',12)
xlabel('multiplying factor')
ylabel('velocity STD mean')
grid on
title(sprintf('LOC-coeff = %4.2e',alpha_coh))

HorizontalSize = 20;
VerticalSize   = 10;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);

set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');


dirprint = '../../slideslastPresentation';
figure(numfig)
printfilename = sprintf('%s/fig1CRBmodif%iLOCfact%iBIS.eps',dirprint,modifloc,not(alpha_coh==0));

fileprintepscmd = ...
    sprintf('print -depsc -loose %s',printfilename);
fileeps2pdf = sprintf('!epstopdf %s',printfilename);
rmeps = sprintf('!rm  %s',printfilename);

if printflag
    eval(fileprintepscmd)
    eval(fileeps2pdf)
    eval(rmeps)
end

figure(numfig2)
printfilename2 = sprintf('%s/fig2CRBmodif%iLOCfact%iBIS.eps',dirprint,modifloc,not(alpha_coh==0));

fileprintepscmd2 = ...
    sprintf('print -depsc -loose %s',printfilename2);

fileeps2pdf2 = sprintf('!epstopdf %s',printfilename2);
rmeps2 = sprintf('!rm  %s',printfilename2);

if printflag
    eval(fileprintepscmd2)
    eval(fileeps2pdf2)
    eval(rmeps2)
end
