%===========================================
function y = delayedsignalF(x,t0)
%===========================================
% delay a signal with a non integer value
% (computation in frequency domain)
%==========
% synopsis:
%     y=delayedsignalF(x,t0)
%==========
% the delay t0 is a non integer delay
% expressed wrt the sampling time:
%   t0 = 1 corresponds to one time dot
%   t0 may be positive, negative, non integer
%   t0>0: shift to the right
%   t0<0: shift to the left
% Rk: the length of FFT is 2^(nextpow2(N)+1
%===========================================

%===============================================
% M. Charbit, Jan. 2010
%===============================================
N         = length(x);
p         = nextpow2(N)+1;
Lfft      = 2^p;
Lffts2    = Lfft/2;
fftx      = fft(x,Lfft);
fftdelay  = ...
    exp(-2j*pi*t0*[(0:Lffts2)'; ...
    -Lfft+((Lffts2+1:Lfft-1)')]/Lfft);
fftdelay(Lffts2+1) = real(fftdelay(Lffts2+1));
ifftdelay          = ifft(fftx .* fftdelay);
y                  = ifftdelay(1:N);
%===============================================

