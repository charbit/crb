clear
allcolors = ['g.';'m.';'r.';'k.';'b.';'rx';'yx';'mx';'rx';'kx';'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

listaz = linspace(0,360,30)';
Laz = length(listaz);
stdaz = zeros(Laz,1);
stdel = zeros(Laz,1);

aec.e_deg      = 45;
aec.c_mps      = 340;
Fs_Hz          = 20;
SNR_dB         = -25;
% T_sec is directly in relationship
% with the max. of delay through the station, 
% which is equal to aec.c_mps * 3000m => 10 sec.
% 10 times this time could be a good choice. 
% But for computation it is almost 
% equivalent to take 10 times less and adjust the
% SNR by the ratio.
T_sec          = 100;
T_cal          = 10;
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
K             = fix(N/2)-1;
frequency_Hz  = (1:K)'*Fs_Hz/N;
sigma2noise   = 10^(-SNR_dB/10);

alpha_coh     = 0.05;
Llistfactor   = 3;
listfactor    = [500, 1200, 5000];linspace(500,2000,Llistfactor);

choix = 3;
switch choix
    case 5
        M               = 3;
        xsensor0        = zeros(M,3);        
        xsensor0(:,2)   = [-1;0;1];
    case 4
        M               = 8;
        xsensor0(:,1:2) = randn(M,2)/3;
        xsensor0(:,3)   = zeros(M,1);        
    case 1
        M               = 8;
        aux             = exp(2j*pi*(0:M-1)'/M);
        xsensor0(:,1)   = real(aux);
        xsensor0(:,2)   = imag(aux);
        xsensor0(:,3)   = 0*ones(M,1);
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
end
    M = size(xsensor0,1);
%     case 1
%         nb_capteurs  = 3;
%         th0          = 0;%pi*rand;
%         P0           = [cos(th0+2*pi*(0:nb_capteurs-1)'/nb_capteurs) ...
%                         sin(th0+2*pi*(0:nb_capteurs-1)'/nb_capteurs)];
%         position_capteurs(:,1:2) = lambda * P0 /sin(pi/nb_capteurs)/4;
%         %P0'*P0/(nb_capteurs/2)
%     case 2
%         position_capteurs = (lambda/2) * [ ...
%             0  1 ; ...
%             0 -1 ; ...
%             -1  0 ];
%     case 3

for ifactor=1:Llistfactor
    factor    = listfactor(ifactor);
    xsensor   = factor * xsensor0;
    % end
    sensordistance = zeros(M,M);
    for im1=1:M
        for im2=1:M
            sensordistance(im1,im2)=norm(xsensor(im1,:)-xsensor(im2,:));
        end
    end
    C     = zeros(K,M,M);
    for ik=1:K
        lambda_k = aec.c_mps/frequency_Hz(ik);
        for im1=1:M
            for im2=1:M
                C(ik,im1,im2)=exp(-alpha_coh*(sensordistance(im1,im2)/lambda_k) );
            end
        end
    end
    for iaz=1:Laz
        aec.a_deg = listaz(iaz);
        
        CRB = evalCRBwithLOC(xsensor, sigma2noise, C, aec, frequency_Hz);
        
        stdaz(iaz) = sqrt(CRB.aec(1,1))*180/pi;
        stdel(iaz) = sqrt(CRB.aec(2,2))*180/pi;
    end
    
    figure(numfig)
    subplot(131)
    if not(ifactor==1)
        hold on
    end
    x = RHO * stdaz .* exp(1j*pi*listaz/180);
    plot(x,'.-','color',allcolors(ifactor))
    hold off
    Mmax = 15;%max(abs(x));
    set(gca,'xlim',Mmax*[-1,1])
    set(gca,'ylim',Mmax*[-1,1])
    axis('square')

    drawnow
    %==
    subplot(132)
    if not(ifactor==1)
        hold on
    end
    x = RHO * stdel .* exp(1j*pi*listaz/180);
    plot(x,'.-','color',allcolors(ifactor))
    hold off
    Mmax = 15;
    set(gca,'xlim',Mmax*[-1,1])
    set(gca,'ylim',Mmax*[-1,1])
    axis('square')
    drawnow
end
subplot(133)
plot(xsensor0(:,1),xsensor0(:,2),'o')
set(gca,'xlim',2*[-1,1])
set(gca,'ylim',2*[-1,1])
axis('square')
% CRB.slowness
% CRB.aec
% sqrt(CRB.aec(1,1))*180/pi
% sqrt(CRB.aec(2,2))*180/pi
% sqrt(CRB.aec(3,3))
