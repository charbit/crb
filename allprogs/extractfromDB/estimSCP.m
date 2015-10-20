function [allSDs, time_sec, frqsFFT_Hz] = ...
    estimSCP(xU,xR,GREF,overlapFFT, ...
    NaverageFFTs, overlapSD, Fs_Hz, smoothwindow)
%==========================================================================
% Perform the spectral components of the two signals xU et xR.
%==========================================================================
% The code uses the Welch's approach. The signal is shared into DFT
% windows of which the length is the length of GREF, with the
% overlap rate of OVERLAPFFT. Then the specral components is averaged
% on NaverageFFTs DFT blocks. Therefore each spectral block corresponds
% to a time period reported in TIME_SEC.
%
% Inputs:
%    xU: signal observed on the SUT (T x 1)
%    xR: signal observed on the SREF (T x 1)
%    GREF: frequency response of the SREF (Lfft x 1)
%    overlapFFT: between 0 and 1, overlap on the FFT block
%    NaverageFFTs: number of FFT-length to averaging
%    overlapSD: between 0 and 1, overlap on the averaging
%                  Spectral Density block
%    Fs_Hz: sampling frequency in Hz
%    smoothwindow: character word as 'rect', 'hann', ...
%            [default: rectangular]
% Outputs:
%    allSDs.RR = auto-spectrum of SREF
%    allSDs.UU = auto-spectrum of SUT
%    allSDs.UR = cross-spectrum of (SUT,SREF)
%    allSDs.MSC = Magnitude Square Coherence
%    allSDs.Rsup = SUU/SUR
%    allSDs.Rinf = SUR/SRR
%    allSDs.det = determinant of the spectral matrix
%
%    time_sec.FFT: grid of time for the FFT blocks (in second)
%    time_sec.SD: grid of time for the Spectral Density blocks (in second)
%    time_sec.signals: grid of time for the samplig period (in second)
%    frqsFFT_Hz: grid of frequency for the interval [0, Fs_Hz] (in Hz)
%==========================================================================
%==========================================================================
xU   = xU(:);
xR   = xR(:);
N    = length(xU);
Lfft = length(GREF);
sqrtLfft    = sqrt(Lfft);
shiftSignal = fix((1-overlapFFT)*Lfft);
NblocksFFT  = fix((N-(Lfft-shiftSignal))/shiftSignal);
allFFTsRR   = zeros(Lfft,NblocksFFT);
allFFTsUU   = zeros(Lfft,NblocksFFT);
if exist('smoothwindow','var')
    switch smoothwindow
        case 'hann'
            Hwin = hann(Lfft,'periodic');
        case 'bartlett'
            Hwin = bartlett(Lfft);
        case 'hamming'
            Hwin = hamming(Lfft);
        case 'rectwin'
            Hwin = rectwin(Lfft);
        case 'blackman'
            Hwin = blackman(Lfft);
    end
else
    Hwin     = hamming(Lfft);
end
%========= normalisation
% not useful if only PSD ratios are considered
Hwin = Hwin *sqrt(Lfft/(Hwin'*Hwin));
for ibF  = 1:NblocksFFT
    ibT  = (ibF-1)*shiftSignal+(1:Lfft);
    xU_i = xU(ibT) .* Hwin;
    xU_i = xU_i-mean(xU_i);
    xR_i = xR(ibT) .* Hwin;
    xR_i = xR_i-mean(xR_i);
    allFFTsUU(:,ibF) = fft(xU_i,Lfft)/sqrtLfft;
    allFFTsRR(:,ibF) = fft(xR_i,Lfft)/sqrtLfft;
end
shiftFFTs = fix((1-overlapSD)*NaverageFFTs);
NSD       = fix(N/Lfft/NaverageFFTs);% fix((NblocksFFT/2-(NaverageFFTs-shiftFFTs))/shiftFFTs);
allSDs    = struct;
time_sec  = struct;
NaverageFFTs_with_overlap = round(NaverageFFTs/(1-overlapFFT)-1);
for ibB=1:NSD,
    indB1 = (ibB-1)*shiftFFTs+1;
    indB2 = indB1+NaverageFFTs_with_overlap-1;
    indB  = fix(indB1):fix(indB2);
    indB  = indB(indB<= NblocksFFT);
    allSDs(ibB).RR  = mean(abs(allFFTsRR(:,indB)) .^ 2,2);
    allSDs(ibB).UU  = mean(abs(allFFTsUU(:,indB)) .^ 2,2);
    allSDs(ibB).UR  = mean(allFFTsRR(:,indB) .* conj(allFFTsUU(:,indB)),2);
    allSDs(ibB).MSC = (abs(allSDs(ibB).UR) .^2) ./...
        (allSDs(ibB).RR .* allSDs(ibB).UU);
    allSDs(ibB).Rsup  = allSDs(ibB).UU ./ allSDs(ibB).UR;
    allSDs(ibB).Rinf  = allSDs(ibB).UR ./ allSDs(ibB).RR;
    allSDs(ibB).det   = (abs(allSDs(ibB).RR) .* abs(allSDs(ibB).UU)) ...
        -(abs(allSDs(ibB).UR) .^2);
end
frqsFFT_Hz = (0:Lfft-1)*Fs_Hz/Lfft;
time_sec.FFT = ((0:NblocksFFT-1)+1/2)*shiftSignal/Fs_Hz;
time_sec.SD =  ((0:NSD-1)+1/2)*shiftFFTs*Lfft/Fs_Hz;
time_sec.signals =  ((0:N-1)+1/2)/Fs_Hz;