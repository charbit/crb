clear
allcolors = ['g.';'m.';'r.';'k.';'b.';'rx';'yx';'mx';'rx';'kx';'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

listaz = linspace(0,360,20)';
Laz = length(listaz);
stdaz = zeros(Laz,1);
stdel = zeros(Laz,1);

aec.e_deg     = 45;
aec.c_mps     = 300;
Fs_Hz         = 20;
sigma2noise   = 1;
K             = 200;
M             = 8;

alpha_coh     = 0.08;
choix = 3;
switch choix
    case 4
        xsensor0(:,1:2) = randn(M,2)/3;
        xsensor0(:,3)   = zeros(M,1);        
    case 1
        aux           = exp(2j*pi*(0:M-1)'/M);
        xsensor0(:,1) = real(aux);
        xsensor0(:,2) = imag(aux);
        xsensor0(:,3) = 0*ones(M,1);
    case 2
        Mon2 = M/2;
        aux1 = exp(2j*pi*(0:Mon2-1)'/Mon2);
        aux2 = exp(2j*pi*(0:Mon2-1)'/Mon2+1j*pi/4)/2;
        xsensor0(:,1) = real([aux1;aux2]);
        xsensor0(:,2) = imag([aux1;aux2]);
        xsensor0(:,3) = 0*ones(M,1);
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
Llistfactor = 3;
listfactor = linspace(300,5000,Llistfactor);

for ifactor=1:Llistfactor
    factor = listfactor(ifactor);
    xsensor   = factor * xsensor0;
    % end
    sensordistance = zeros(M,M);
    for im1=1:M
        for im2=1:M
            sensordistance(im1,im2)=norm(xsensor(im1,:)-xsensor(im2,:));
        end
    end
    fK    = (1:K)'*Fs_Hz/2/K;
    C     = zeros(K,M,M);
    for ik=1:K
        lambda_k = aec.c_mps/fK(ik);
        for im1=1:M
            for im2=1:M
                C(ik,im1,im2)=exp(-alpha_coh*(sensordistance(im1,im2)/lambda_k) );
            end
        end
    end
    
    for iaz=1:Laz
        aec.a_deg = listaz(iaz);
        
        CRB = CRBcoherence(xsensor, sigma2noise, C, aec, K, Fs_Hz);
        
        stdaz(iaz) = sqrt(CRB.aec(1,1))*180/pi;
        stdel(iaz) = sqrt(CRB.aec(2,2))*180/pi;
    end
    
    figure(1)
    subplot(121)
    if not(ifactor==1)
        hold on
    end
    x = stdaz .* exp(1j*pi*listaz/180);
    plot(x,'.-','color',allcolors(ifactor))
    hold off
    drawnow
    %==
    subplot(122)
    if not(ifactor==1)
        hold on
    end
    x = stdel .* exp(1j*pi*listaz/180);
    plot(x,'.-','color',allcolors(ifactor))
    hold off
    drawnow
end
% CRB.slowness
% CRB.aec
% sqrt(CRB.aec(1,1))*180/pi
% sqrt(CRB.aec(2,2))*180/pi
% sqrt(CRB.aec(3,3))
