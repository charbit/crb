clear
aec.a_deg = 30;
aec.e_deg = 20;
aec.c_mps = 330;

cosa  = cos(aec.a_deg*pi/180);
sina  = sin(aec.a_deg*pi/180);
cose  = cos(aec.e_deg*pi/180);
sine  = sin(aec.e_deg*pi/180);
c_mps = aec.c_mps;
v_mps = c_mps/cose;

slowness_spm    = zeros(3,1);
slowness_spm(1) = cosa*cose/c_mps;
slowness_spm(2) = sina*cose/c_mps;
slowness_spm(3) = sine/c_mps;


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

factor          = 1000;
xsensor         = factor * xsensor0;
delay_s         = xsensor * slowness_spm;

Fs_Hz           = 20;
SNR_dB          = -5;
T_sec           = 100;
N               = fix(T_sec*Fs_Hz);
sigmanoise     = 10^(-SNR_dB/20);

st0             = randn(N,1);
stdelayed       = zeros(N,M);
xt              = zeros(N,M);
for im=1:M
    stdelayed(:,im) = delayedsignalF(st0,delay_s(im)*Fs_Hz);
end

Lruns           = 300;
azimuthhat_deg  = zeros(Lruns,1);
velocityhat_mps = zeros(Lruns,1);
for irun=1:Lruns
    nt              = sigmanoise*randn(N,M);
    xt              = stdelayed+nt;
    
    [xalign, tkl_pts] = alignmentwrt1(xt,1,N);
    tkl_sec = (tkl_pts/Fs_Hz);
    slownesshat_spm = pinv(xsensor)*tkl_sec;
    azimuthhat_deg(irun) = ...
        atan2(slownesshat_spm(2),slownesshat_spm(1))*180/pi;
    velocityhat_mps(irun) = 1 / norm(slownesshat_spm(1:2));
    % [azimuthhat_deg, aec.a_deg]
    % [v_mps    velocityhat_mps]
end
%%
figure(2)
plot(azimuthhat_deg,velocityhat_mps,'.')
set(gca,'xlim',aec.a_deg+1*[-1,1])
set(gca,'ylim',v_mps+1*[-1,1])

