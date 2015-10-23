clear

% The spectral matrix model depends on 6 parameters
% of the covariance matrix. The log-MSC is linear 
% wrt the 6 parameters. Therefore we use least square
% approach.
%
%
% used functions
%   - spectralmatrixestimation
%   - estimSigmatheta
%
%
directorydatafromIDC  = '../AAdataI37/';
filenames              = dir(sprintf('%ssta*240.mat',directorydatafromIDC));

load('../sensorlocation/I37.mat');
xsensors_km   = xsensors_m.coordinates/1000;
Msensors            = size(xsensors_km,1);
combi        = Msensors*(Msensors-1)/2;
cp = 0;
distances = zeros(combi,1);
for i1=1:Msensors-1
    for i2=i1+1:Msensors
        cp=cp+1;
        distances(cp) = norm(xsensors_m.coordinates(i1,:)-...
            xsensors_m.coordinates(i2,:));
    end
end
[sortdistances, indsortdistance] = sort(distances);


%==== filter
Fs_Hz                  = 20;
Norderfilter           = 4;

freqrange_Hz           = [0.06 0.14];

[forward,  reverse]    = butter(Norderfilter,2*freqrange_Hz/Fs_Hz);

signals = cell(Msensors,1);
stime = zeros(Msensors,1);
etime = zeros(Msensors,1);

for im=1:Msensors
    filename1_ii = filenames(im).name;
    cdload = sprintf('sig = load(''%s%s'');',directorydatafromIDC,filename1_ii);
    eval(cdload)
    Lrecords   = length(sig.records);
    LL_max=0;
    for ir =1:Lrecords
        if length(sig.records{ir}.data)>LL_max
            irsave = ir;
            LL_max=length(sig.records{ir}.data);
        end
    end
    signals{im} = [sig.records{irsave}.data];
    stime(im) = sig.records{irsave}.stime;
    etime(im) = sig.records{irsave}.etime;
end
stimeMAX = max(stime);
etimeMIN = min(etime);
signalsproc = cell(Msensors,1);

for im=1:Msensors
    Lim = length(signals{im});
    if stime(im)<=stimeMAX, ds = fix((stimeMAX-stime(im))*Fs_Hz)+1; end
    if etime(im)>=etimeMIN, de = fix((etime(im)-etimeMIN)*Fs_Hz); end
    signalsproc{im} = signals{im}(ds:Lim-de);
end
LSIG =size([signalsproc{:}],1);
SIG = zeros(LSIG,Msensors);
for im=1:Msensors
    SIG(:,im) = [signalsproc{im}];
end
clear signalsproc
clear signals
sigfil                 = SIG;
for im = 1:Msensors
    sigfil(:,im)       = filtfilt(forward,reverse,SIG(:,im));
end

cp=0;
allMSC = cell(combi,1);
for im1 = 1:Msensors-1
    for im2 = im1+1:Msensors
        cp=cp+1;
        taustationary_sec = 500 ;
        ratioDFT2SCP = 5;
        Tfft_sec     = taustationary_sec/ratioDFT2SCP;
        Lfft         = fix(Tfft_sec*Fs_Hz);
        GREF         = ones(Lfft,1);
        [allSDs, time_sec, frqsFFT_Hz] = ...
            estimSDs(sigfil(:,im1),sigfil(:,im2),GREF,0.5, ...
            ratioDFT2SCP, 0, Fs_Hz, 'hann');
        allMSC{cp} = [allSDs.MSC];
    end
end
allMSCsort = allMSC(indsortdistance);
%%
freqrange_pts = fix(freqrange_Hz*Lfft/Fs_Hz);
fqnormrange   = freqrange_pts(1):freqrange_pts(2);
Nfreq         = freqrange_pts(2)-freqrange_pts(1)+1;
Nscp = size(allMSCsort{1},2);
atime = fix(rand*Nscp)+1;
MSCcp = zeros(combi,Nfreq);
for ip=1:combi
    MSCcp(ip,:) = allMSC{ip}(fqnormrange,atime);
end

pi2           = pi*pi;
dSigmaTheta2  = 3;
%=================
% we may verify that the time alignment does not
% change the estimation of SigmaTheta2
% xalign = alignmentwrt1(xn,1,N,Fs_Hz);
% xn = xalign.sig;
%=================
% estimation of the model parameter
[Sigmatheta2, logSkcp_pred, logresidue] = ...
    estimSigmatheta(MSCcp, xsensors_km, fqnormrange, Fs_Hz, 2);

% calculation of the gaussian model
logcoh     = zeros(combi,Nfreq);
freq_Hz    = Fs_Hz*fqnormrange;
freq2_Hz   = freq_Hz .^2;
cp         = 0;
for im1 = 1:Msensors-1
    for im2 = im1+1:Msensors
        cp        = cp+1;
        diffloc_km = xsensors_km(im2,1:2)-xsensors_km(im1,1:2);
        logcoh(cp,:) = -2*pi2*freq2_Hz*(diffloc_km*Sigmatheta2*diffloc_km');
    end
end

% [sortdistances, indsortdistance]
figure(1)
Nselected = 5;

plot(sortdistances/lambda_km/1000, mean(exp(logcoh(:,:)),2),'x-r')
hold on
plot(sortdistances/lambda_km/1000, mean(MSCcp(:,:),2),...
    'ob','markerf','b')
hold off
set(gca,'ylim',[0 1])
xlabel('DWR')
ylabel('MSC')
grid on
% txt=sprintf('T = %i s, Fs = %i Hz\nAnalysis freq. range = [%4.1f %4.1f] Hz\nSelected freq. = %4.1f Hz',...
%     T_sec, Fs_Hz, freqrange_Hz(1), freqrange_Hz(2), selectedfreq_Hz);
% title(txt)

hfig=20;
vfig=13;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
% set(gcf,'position',[8 10 hfig vfig]);
set(gcf,'paperposition',[0 0 hfig vfig]);
% print -depsc /Users/charbit/maurice/etudes/InfrasonCEA/contrat2012/livrable2013/slides/analyseCoherIS31
