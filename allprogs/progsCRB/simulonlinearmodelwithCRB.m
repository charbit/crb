%====================== simulonlinearmodelwithCRB.m
clear
% this program simulates signals with full coherence thanks to
% only pure delays. Then by sub-optimal method which consists 
% of a time alignment followed by pseudo-inverse the sensor 
% location matrix, we estimate the variance of the azimuth
%======
% We also provides the Cramer-Rao bound which is in
% good agreement with the estimated value.
%
% We mainly consider 'flat/flat station' station 
% In this case we consider horizontal velocity
%
%
%

addpath toolbox
aec.a_deg = 10;
aec.e_deg = 30;
aec.c_mps = 340;

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

Fs_Hz           = 20;
SNR_dB          = -5;
T_sec           = 50;
N               = fix(T_sec*Fs_Hz);
K               = fix(N/2);
sigmanoise      = 10^(-SNR_dB/20);
%==
stationselect = 6;
switch stationselect
    case 5
        M               = 3;
        xsensor0        = zeros(M,3);
        xsensor0(:,1)   = [-1;0;1];
    case 4
        M               = 4;
        xsensor0(:,1:2) = randn(M,2)/3;
        xsensor0(:,3)   = zeros(M,1);
    case 1
        M               = 8;
        aux             = exp(2j*pi*(0:M-1)'/M);
        xsensor0(:,1)   = real(aux);
        xsensor0(:,2)   = imag(aux);
        xsensor0(:,3)   = 0;
    case 2
        M               = 8;
        Mon2            = M/2;
        aux1            = exp(2j*pi*(0:Mon2-1)'/Mon2);
        aux2            = exp(2j*pi*(0:Mon2-1)'/Mon2+1j*pi/4)/2;
        xsensor0(:,1)   = real([aux1;aux2]);
        xsensor0(:,2)   = imag([aux1;aux2]);
        xsensor0(:,3)   = 0;
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
    case 6
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
end
M               = size(xsensor0,1);
factor          = 1000;
xsensor         = factor * xsensor0;
delay_s         = xsensor * slowness_spm;

st0             = randn(N,1);
stdelayed       = zeros(N,M);
xt              = zeros(N,M);
for im=1:M
    stdelayed(:,im) = delayedsignalF(st0,delay_s(im)*Fs_Hz);
end
coherencematrix = NaN;
CRB = evalCRBwithLOC(xsensor,sigmanoise^2,...
    coherencematrix,aec,Fs_Hz*(1:K)'/N);

%===== simulation =====
Lruns           = 500;
azimuthhat_deg  = zeros(Lruns,1);
velocityhat_mps = zeros(Lruns,1);
for irun=1:Lruns
    nt                    = sigmanoise*randn(N,M);
    xt                    = stdelayed+nt;   
    [xalign, tkl_pts]     = alignmentwrt1(xt,1,N,4);
    tkl_sec               = (tkl_pts/Fs_Hz);
    slownesshat_spm       = pinv(xsensor)*tkl_sec;
    azimuthhat_deg(irun)  = ...
        atan2(slownesshat_spm(2),slownesshat_spm(1))*180/pi;
    velocityhat_mps(irun) = 1 / norm(slownesshat_spm(1:2));
end
hatnoiseratio = sqrt(CRB.av(1,1))*180/pi/std(azimuthhat_deg);
% sqrt(CRB.av(2,2))/ std(velocityhat_mps)


%%
figure(4)
plot(azimuthhat_deg,velocityhat_mps,'.')
set(gca,'xlim',aec.a_deg+0.3*[-1,1])
set(gca,'ylim',v_mps+2*[-1,1])
hold on
X0 = [aec.a_deg; v_mps];
R = CRB.av;
R(1,1) = CRB.av(1,1)*180*180/pi/pi;
R(1,2) = CRB.av(1,2)*180/pi;
R(2,1) = R(1,2);
R(2,2) = CRB.av(2,2);
alphapercent = 0.95;
ellipse(X0, R/hatnoiseratio/hatnoiseratio, alphapercent, 'r');
hold off
