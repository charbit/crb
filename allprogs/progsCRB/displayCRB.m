clear
allcolors = ['g.';'m.';'r.';'k.';'b.';'rx';'yx';'mx';'rx';'kx';'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

listaz = linspace(0,360,50)';
Laz = length(listaz);
stdaz = zeros(Laz,1);
stdel = zeros(Laz,1);
stdvel = zeros(Laz,1);

aec.e_deg      = 20;
aec.c_mps      = 340;
Fs_Hz          = 20;
SNR_dB         = 20;
% T_sec is directly in relationship
% with the max. of delay through the station,
% which is equal to aec.c_mps * 3000m => 10 sec.
% 10 times this time could be a good choice.
% But for computation it is almost
% equivalent to take 10 times less and adjust the
% SNR by the ratio.
T_sec          = 100;
T_cal          = 100;
% Then we correct by sqrt(T_cal/T_sec)
RHO            = sqrt(T_cal/T_sec);
N              = fix(T_cal*Fs_Hz);
%==================================
% Rk: to use Whittle's formula
% the frequencies must be those of the DFT for
% N signal samples. Except 0 and N/2 if N/2 is
% integer, because the DFT is real. If the signal
% is real, the DFT has hermitian symmetry, therefore
% we have to take only the half part of the
% frequency domain
% If we use a filter bank, we have to decimate and
% perform the DFT in each band. For example
% if we have a band between 0.1Hz and 0.8Hz, i.e.
% a bandwidth of 0.7 Hz, we have to re-sample the
% 20 Hz signals at 1.4Hz. Instead of N=2000 points
% we have now N=140 points.
%
K             = 18;%fix(N/2)-1;
frequency_Hz  = (1:K)'*Fs_Hz/N;
sigma2noise   = 10^(-SNR_dB/10);

alpha_coh     = 8e-5;%8e-5;%0;%1.42e-4;%0.008;
Llistfactor   = 12;
listfactor    = linspace(500,2000,Llistfactor);
% listfactor    = [500, 800, 1000, 1500, 2000];
% Llistfactor   = length(listfactor);

choice = 6;
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
        %
end
numfig         = 1;

if 0
    xsensors_new_m = transform2isotrop(xsensors_centered_norm_m,1.2);
else
    xsensors_new_m = xsensors_centered_norm_m;
end
M = size(xsensors_new_m,1);
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
            sensordistance(im1,im2)=norm(xsensors_fact_m(im1,:)-xsensors_fact_m(im2,:));
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
        CRB          = evalCRBwithLOC(xsensors_fact_m, sigma2noise, C, polparameters, frequency_Hz);
        stdaz(iaz)   = sqrt(CRB.av(1,1))*180/pi;
        stdvel(iaz)  = sqrt(CRB.av(2,2));
        
    end
    azasaz = RHO * stdaz .* exp(1j*pi*listaz/180);
    meanstdaz(ifactor) = mean(stdaz);
            velasaz = RHO * stdvel .* exp(1j*pi*listaz/180);

        meanstdvel(ifactor) = mean(stdvel);

        figure(numfig)
        subplot(131)
    % Mmax = 2;
    % set(gca,'xlim',Mmax*[-1,1])
    % set(gca,'ylim',Mmax*[-1,1])
    % set(gca,'xtick',[0],'ytick',[0])
    xlabel('azimuth');
    axis('square')
    grid on
    drawnow
    
    subplot(132)
    % Mmax = 10;
    % set(gca,'xlim',Mmax*[-1,1])
    % set(gca,'ylim',Mmax*[-1,1])
    % set(gca,'xtick',[0],'ytick',[0])
    xlabel('hor. velocity');
    axis('square')
    grid on
    drawnow
    
    subplot(133)
    plot(xsensors_fact_m(:,1),xsensors_fact_m(:,2),'o')
    subplot(133)
    % Mmax = 1;
    % set(gca,'xlim',Mmax*[-1,1])
    % set(gca,'ylim',Mmax*[-1,1])
    axis('square')
    % set(gca,'xtick',[0],'ytick',[0])
    grid on
    % CRB.slowness
    % CRB.aec
    % sqrt(CRB.aec(1,1))*180/pi
    % sqrt(CRB.aec(2,2))*180/pi
    % sqrt(CRB.aec(3,3))
    
    subplot(131)
        if not(ifactor==1)
            hold on
        end
        plot(azasaz,'.-','color',allcolors(ifactor))
        hold off
    
        %==
        subplot(132)
        if not(ifactor==1)
            hold on
        end
        plot(velasaz,'.-','color',allcolors(ifactor))
        hold off
end

subplot(131)
% Mmax = 2;
% set(gca,'xlim',Mmax*[-1,1])
% set(gca,'ylim',Mmax*[-1,1])
% set(gca,'xtick',[0],'ytick',[0])
xlabel('azimuth');
axis('square')
grid on
drawnow

subplot(132)
% Mmax = 10;
% set(gca,'xlim',Mmax*[-1,1])
% set(gca,'ylim',Mmax*[-1,1])
% set(gca,'xtick',[0],'ytick',[0])
xlabel('hor. velocity');
axis('square')
grid on
drawnow

subplot(133)
plot(xsensors_fact_m(:,1),xsensors_fact_m(:,2),'o')
subplot(133)
% Mmax = 1;
% set(gca,'xlim',Mmax*[-1,1])
% set(gca,'ylim',Mmax*[-1,1])
axis('square')
% set(gca,'xtick',[0],'ytick',[0])
grid on
% CRB.slowness
% CRB.aec
% sqrt(CRB.aec(1,1))*180/pi
% sqrt(CRB.aec(2,2))*180/pi
% sqrt(CRB.aec(3,3))
%
figure(numfig)
HorizontalSize = 24;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
% set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);

set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
set(gcf, 'InvertHardCopy', 'off');
subplot(132)
% title(sprintf('LOC-coeff = %4.2e\n%s',listlegs))

fileprintepscmd = sprintf('print -depsc -loose ../../figures/CRBI0%i.eps',stationnumber);
% eval(fileprintepscmd)



%%
numfig2 = 3;
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

set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
set(gcf, 'InvertHardCopy', 'off');

fileprintepscmd = sprintf('print -depsc -loose ../../figures/CRBI0%i.eps',stationnumber);
% eval(fileprintepscmd)

