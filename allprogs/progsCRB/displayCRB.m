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
SNR_dB         = 0;
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
numfig         = 1;
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

alpha_coh     = 1.42e-4;%0.008;
Llistfactor   = 3;
listfactor    = [500, 1000, 2000];%linspace(500,2000,Llistfactor);
Llistfactor   = length(listfactor);

choice = 6;
switch choice
    case 5
        M               = 8;
        xsensor0        = zeros(M,3);
        xsensor0(:,1)   = (-M+1:2:M-1)'/M;
    case 4
        M               = 8;
        xsensor0(:,1:2) = randn(M,2)/3;
        xsensor0(:,3)   = zeros(M,1);
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
    case 3
        xsensor0 = [ ...
            -0.05580864939648136, 0.1876414387122062,   0; ...
            0.2276638554764773,  0.08756600473273159,  0; ...
            0.1213616661490549, -0.2126602972080624,   0; ...
            -0.1266767756150987, -0.09034587789893783,  0; ...
            -0.02746139890923584, 0.009729556081326833, 0; ...
            0.922171492415503,   0.7102575939423894,   0; ...
            0.1851429797456091, -1.068861232371144,    0; ...
            -1.246393169863814,   0.3766728140087006,   0];
    case 6  %I37
        xsensor0 = 0.1*[ ...
            0                   0                   0; ...
            -0.046877141457517   1.955947627983987  -0.003000000000000; ...
            2.294699908569455  -0.292623042389750   0.004000000000000; ...
            -1.529189879004843  -1.626142297079787  -0.012000000000000; ...
            4.450550248185173   5.573021452268586   0.044000000000000; ...
            0.396748844208196   9.435171689242125  -0.006000000000000; ...
            9.679406059307512   5.422483168104663   0.061000000000000; ...
            3.943046386367641  -5.037650076216087  -0.118000000000000; ...
            -6.480990306646564  -5.857643549628555  -0.017000000000000; ...
            -6.833279994535260   1.702617467129603  -0.059000000000000];
        
        %
end

if 0
    gX          = xsensor0';
    xsensor0new = transform2isotrop(xsensor0);
else
    xsensor0new = xsensor0;
end

M = size(xsensor0new,1);
for ifactor=1:Llistfactor
    factor    = listfactor(ifactor);
    xsensor_m   = factor * xsensor0new;
    % end
    sensordistance = zeros(M,M);
    for im1=1:M
        for im2=1:M
            sensordistance(im1,im2)=norm(xsensor_m(im1,:)-xsensor_m(im2,:));
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
        
        CRB = evalCRBwithLOC(xsensor_m, sigma2noise, C, aec, frequency_Hz);
        
        stdaz(iaz)  = sqrt(CRB.av(1,1))*180/pi;
        stdvel(iaz)  = sqrt(CRB.av(2,2));
    end
    
    figure(numfig)
    subplot(131)
    if not(ifactor==1)
        hold on
    end
    x = RHO * stdaz .* exp(1j*pi*listaz/180);
    plot(x,'.-','color',allcolors(ifactor))
    hold off    
    
    %==
    subplot(132)
    if not(ifactor==1)
        hold on
    end
    x = RHO * stdvel .* exp(1j*pi*listaz/180);
    plot(x,'.-','color',allcolors(ifactor))
    hold off
end
subplot(131)
% Mmax = 1.8;
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
set(gca,'xtick',[0],'ytick',[0])
xlabel('hor. velocity');
axis('square')
grid on
drawnow



subplot(133)
plot(xsensor0new(:,1),xsensor0new(:,2),'o')
subplot(133)
Mmax = 1.5;
set(gca,'xlim',Mmax*[-1,1])
set(gca,'ylim',Mmax*[-1,1])
axis('square')
set(gca,'xtick',[0],'ytick',[0])
grid on
% CRB.slowness
% CRB.aec
% sqrt(CRB.aec(1,1))*180/pi
% sqrt(CRB.aec(2,2))*180/pi
% sqrt(CRB.aec(3,3))
%%
figure(1)
HorizontalSize = 24;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);

set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
set(gcf, 'InvertHardCopy', 'off');
subplot(132)
title(sprintf('LOC-coeff = %4.2e\ngreen: %ix, magenta: %ix, red: %ix',alpha_coh,listfactor))

stationnumber=37;
fileprintepscmd = sprintf('print -depsc -loose ../../figures/CRBI0%i.eps',stationnumber);
% eval(fileprintepscmd)