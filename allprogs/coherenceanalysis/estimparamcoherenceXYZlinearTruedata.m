clear
switch computer
    case 'GLNXA64'
        addpath /dvlscratch/SHI/users/charbit/ProjectIMS2015b/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
otherwise
    addpath /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
end

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
filenames              = dir(sprintf('%s*.mat',directorydatafromIDC));

load('../sensorlocation/I37.mat');
xsensors_m   = xsensors_m.coordinates;
Msensors            = size(xsensors_m,1);
combi        = Msensors*(Msensors-1)/2;
cp = 0;
distances = zeros(combi,1);
for i1=1:Msensors-1
    for i2=i1+1:Msensors
        cp=cp+1;
        distances(cp) = norm...
            (xsensors_m(i1,:)-xsensors_m(i2,:));
    end
end
[sortdistances, indsortdistance] = sort(distances);

cp = 0;
directions = zeros(combi,1);
for i1=1:Msensors-1
    for i2=i1+1:Msensors
        cp=cp+1;
        directions(cp) = (180/pi)*atan(...
            (xsensors_m(i1,2)-xsensors_m(i2,2)) / ...
            (xsensors_m(i1,1)-xsensors_m(i2,1)));
    end
end
[sortdirections, indsortdirection] = sort(directions);


%==== filter
Fs_Hz                  = 20;
Norderfilter           = 2;

freqrange_Hz           = [0.08 0.12];

[forward,  reverse]    = butter(Norderfilter,2*freqrange_Hz/Fs_Hz);


selectednumfile = 3;

filename1_ii = filenames(selectednumfile).name;
cdload = sprintf('load(''%s%s'');',directorydatafromIDC,filename1_ii);
eval(cdload)


sigfil                 = observations.data;
for im = 1:Msensors
    SIG = observations.data(:,im);
    sigfil(:,im)       = filtfilt(forward,reverse,SIG);
end
T_pts = size(sigfil,1);
T_sec = T_pts/Fs_Hz;
cp=0;
allMSC = cell(combi,1);
for im1 = 1:Msensors-1
    for im2 = im1+1:Msensors
        cp=cp+1;
        taustationary_sec = 1000 ;
        ratioDFT2SCP      = 10;
        Tfft_sec          = taustationary_sec/ratioDFT2SCP;
        Lfft              = fix(Tfft_sec*Fs_Hz);
        GREF              = ones(Lfft,1);
        [allSDs, time_sec, frqsFFT_Hz] = ...
            estimSCP(sigfil(:,im1),sigfil(:,im2),GREF,0.5, ...
            ratioDFT2SCP, 0, Fs_Hz, 'hann');
        allMSC{cp} = [allSDs.MSC];
    end
end
%%
% %%
% figure(1)
% clf
% for iss=1:45,
%     subplot(9,5,iss),
%     pcolor(time_sec.SD/3600,frqsFFT_Hz, allMSCsort{iss}),
%     shading flat,
%     set(gca,'xticklabel',[])
%      set(gca,'ylim',[0, 0.14])
%      if iss>1
%     set(gca,'yticklabel',[])
%      end
% end
% 


atime=0;
%%
% 
atime = atime+1;

freqrange_bins = fix(freqrange_Hz*Lfft/Fs_Hz);
indfqrange     = freqrange_bins(1):freqrange_bins(2);
Nfreq          = freqrange_bins(2)-freqrange_bins(1)+1;
Nscp           = size(allMSC{1},2);
% selectdirections=find(or(and(directions>-90,directions<-70), ...
%     and(directions>70,directions<90)));
% Lselectdirections = length(selectdirections)
MSCinBW = zeros(combi,Nfreq);
for ip=1:combi
    MSCinBW(ip,:) = allMSC{ip}(indfqrange,atime);
end
pi2           = pi*pi;
%=================
% we may verify that the time alignment does not
% change the estimation of SigmaTheta2
% xalign = alignmentwrt1(xn,1,N,Fs_Hz);
% xn = xalign.sig;
%=================
% MSCinBW = 0.3*rand(size(MSCinBW));
% estimation of the model parameter
[Sigmatheta2_s2pm2, logSkcp_pred, logresidue] = ...
    estimSigmatheta(MSCinBW, xsensors_m, ...
    frqsFFT_Hz(indfqrange)/Fs_Hz, Fs_Hz, 2,0.000000);

% calculation of the gaussian model
predictMSC     = zeros(combi,Nfreq);
freq_Hz    = frqsFFT_Hz(indfqrange);
freq2_Hz   = freq_Hz .^2;
cp         = 0;
for im1 = 1:Msensors-1
    for im2 = im1+1:Msensors
        cp        = cp+1;
        diffloc_m = xsensors_m(im2,1:2)-xsensors_m(im1,1:2);
        predictMSC(cp,:) = ...
            exp(-4*pi2*freq2_Hz*(diffloc_m*Sigmatheta2_s2pm2*diffloc_m'));
    end
end

% [sortdistances, indsortdistance]
figure(1)

sortpredictMSC = predictMSC(indsortdistance,:);
sortMSCinBW    = MSCinBW(indsortdistance,:);

plot(sortdistances,  median(sortpredictMSC,2),'x-r')
hold on
plot(sortdistances, median(sortMSCinBW,2),'o-b','markerf','b')
plot(sortdistances, sortMSCinBW,'.','color',0.7*ones(3,1))
hold off
% set(gca,'ylim',[0 1])
xlabel('DWR')
ylabel('MSC')
grid on
txt=sprintf('T = %4.1f H, Fs = %i Hz\nAnalysis freq. range = [%4.2f %4.2f] Hz',...
    T_sec/3600, Fs_Hz, freqrange_Hz(1), freqrange_Hz(2));
title(txt)

hfig=20;
vfig=13;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
% set(gcf,'position',[8 10 hfig vfig]);
set(gcf,'paperposition',[0 0 hfig vfig]);
% print -depsc /Users/charbit/maurice/etudes/InfrasonCEA/contrat2012/livrable2013/slides/analyseCoherIS31
