%=========== displayCRBlin.m
% this program computes the CRB for the linear model
% associated to the M-1 delays w.r.t. to the sensor 1 and 
% the difference of sensor locations.
% Usually these delays are performed by cross-correlation
% maximization. We build the M-1 x d matrix H of the 
% difference of locations wrt the sensor 1 and computes
% inv(H'*H). To transform in a,e,c we use the Jacobian.
%
clear
listaz = linspace(0,360,100)';
Laz = length(listaz);
stdaz = zeros(Laz,1);
stdel = zeros(Laz,1);
stdvel = zeros(Laz,1);

transform2isoflag = 0;
aec.e_deg      = 20;
aec.c_mps      = 340;

choix = 6;
switch choix
    case 5
        M               = 8;
        xsensors0_m        = zeros(M,3);
        xsensors0_m(:,1)   = (-M+1:2:M-1)'/M;
        xsensors0_m(:,2)   = xsensors0_m(:,1);
    case 4
        M               = 8;
        xsensors0_m(:,1:2) = randn(M,2)/3;
        xsensors0_m(:,3)   = zeros(M,1);
    case 1
        M               = 4;
        aux             = exp(2j*pi*(0:M-1)'/M);
        xsensors0_m(:,1)   = real(aux);
        xsensors0_m(:,2)   = imag(aux);
        xsensors0_m(:,3)   = 1*ones(M,1);
    case 2
        M               = 8;
        Mon2            = M/2;
        aux1            = exp(2j*pi*(0:Mon2-1)'/Mon2);
        aux2            = exp(2j*pi*(0:Mon2-1)'/Mon2+1j*pi/4)/2;
        xsensors0_m(:,1)   = real([aux1;aux2]);
        xsensors0_m(:,2)   = imag([aux1;aux2]);
        xsensors0_m(:,3)   = 0*ones(M,1);
    case 3
        xsensors0_m = [ ...
            -0.05580864939648136, 0.1876414387122062,   0; ...
            0.2276638554764773,  0.08756600473273159,  0; ...
            0.1213616661490549, -0.2126602972080624,   0; ...
            -0.1266767756150987, -0.09034587789893783,  0; ...
            -0.02746139890923584, 0.009729556081326833, 0; ...
            0.922171492415503,   0.7102575939423894,   0; ...
            0.1851429797456091, -1.068861232371144,    0; ...
            -1.246393169863814,   0.3766728140087006,   0];
    case 6
        xsensors0_m = 0.1*[ ...
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
        
end
% if 1 we transform the array in isotrope array
% for coherent signal.
if transform2isoflag
    newX     = transform2isotrop(xsensors0_m);
    xsensors0_mnew = newX;
else
    xsensors0_mnew = xsensors0_m;
end
% 
M = size(xsensors0_mnew,1);

H3 = zeros(M-1,3);
for im=2:M
    H3(im-1,:) = (xsensors0_mnew(im,:)-xsensors0_mnew(1,:));
end
H3H3t = (H3'*H3);

for iaz=1:Laz
    aec.a_deg = listaz(iaz);
    
    
    cosa = cos(aec.a_deg*pi/180);
    sina = sin(aec.a_deg*pi/180);
    cose = cos(aec.e_deg*pi/180);
    sine = sin(aec.e_deg*pi/180);
    c_mps    = aec.c_mps;
    v_mps = c_mps/cose;
    
    Jacobaec_k = ([...
        -sina*cose/c_mps -cosa*sine/c_mps -cosa*cose/c_mps/c_mps; ...
        cosa*cose/c_mps -sina*sine/c_mps -sina*cose/c_mps/c_mps;...
        0 cose/c_mps -sine/c_mps/c_mps]);
    
    Jacobav_k = ([...
        -sina/v_mps -cosa/v_mps/v_mps 0; ...
        cosa/v_mps -sina/v_mps/v_mps 0;...
        0 -sine/cose/v_mps/v_mps 1/cose/cose/v_mps ...
        ]);
        
    CRB.linaec = (Jacobaec_k\H3H3t)/Jacobaec_k';
    CRB.linav  = (Jacobav_k\H3H3t)/Jacobav_k';
    
    figure(1)
    subplot(222)
    ellipse(zeros(2,1),CRB.linav(1:2,1:2),0.95);
    hold on

    stdaz(iaz)   = sqrt(CRB.linav(1,1))*180/pi;
    stdvel(iaz)  = sqrt(CRB.linav(2,2));
    
end
figure(1)
subplot(222)
axis('square')
set(gca,'xtick',[],'ytick',[])
hold off
title('BCR on (a,v)')

subplot(223)
x = stdaz .* exp(1j*pi*listaz/180) / max(stdaz);
plot(x,'.-')
hold off
set(gca,'xlim',[-1,1])
set(gca,'ylim',[-1,1])
set(gca,'xtick',[],'ytick',[])
axis('square')
title('incertitude in azimut')

subplot(224)
x = stdvel .* exp(1j*pi*listaz/180) / max(stdvel);
plot(x,'.-')
hold off
set(gca,'xlim',[-1,1])
set(gca,'ylim',[-1,1])
set(gca,'xtick',[],'ytick',[])
axis('square')
title('incertitude in trace velocity')

subplot(221)
plot(xsensors0_m(:,1),xsensors0_m(:,2),'o')
set(gca,'xlim',2*[-1,1])
set(gca,'ylim',2*[-1,1])
set(gca,'xtick',[],'ytick',[])
axis('square')
title('sensor locations')

