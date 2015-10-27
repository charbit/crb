%====================== estimationwithFB.m =============================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];
switch computer
    case 'GLNXA64'
        addpath /dvlscratch/SHI/users/charbit/ProjectIMS2015b/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
otherwise
    addpath /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
end

stationnumber         = 37;
directorydatafromIDC  = sprintf('../../../AAdataI%i/',stationnumber);

filenames              = dir(sprintf('%s*.mat',directorydatafromIDC));
nbfiles                = length(filenames);
Fs_Hz                  = 20;
Ts_sec                 = 1/Fs_Hz;

ifile = 5;
filename1_ii = filenames(ifile).name;
cdload = sprintf('load(''%s%s'');',directorydatafromIDC,filename1_ii);
eval(cdload)

xsensors_m = observations.xsensors_m.coordinates;
Msensors   = size(xsensors_m,1);
cp = 0;
combi = Msensors*(Msensors-1)/2;
distances = zeros(combi,1);
for i1=1:Msensors-1
    for i2=i1+1:Msensors
        cp=cp+1;
        distances(cp) = norm(xsensors_m(i1,:)-xsensors_m(i2,:));
    end
end
[sortdistances, indsortdistance] = sort(distances);
% orientations = zeros(combi,1);
% for i1=1:Msensors-1
%     for i2=i1+1:Msensors
%         cp=cp+1;
%         orientations(cp) = atan((xsensors_m(i1,2)-xsensors_m(i2,2)) / >>>
%         (xsensors_m(i1,2)-xsensors_m(i2,2));
%     end
% end
% [sortdistances, indsortdistance] = sort(distances);

signals = cell(Msensors,1);
stime = zeros(Msensors,1);
etime = zeros(Msensors,1);


LSIG = size(observations.data,1);
SIG = zeros(LSIG,Msensors);
for im=1:Msensors
    SIG(:,im) = observations.data(:,im);
end

allMSC = cell(combi,1);

bandwidthFilter_Hz  = [0.02 0.2];
[forward,  reverse] = butter(2,2*bandwidthFilter_Hz/Fs_Hz);
% for im = 1:Msensors
%     SIG(:,im)       = filtfilt(forward,reverse,SIG(:,im));
% end
bandwidthMSC_Hz   = [0.06 0.16];
taustationary_sec = 500 ;
ratioDFT2SCP      = 5;

%=====================
%=====================
%=====================
cp=0;
for im1 = 1:Msensors-1
    for im2 = im1+1:Msensors
        cp = cp+1;
        signal1_centered=SIG(:,im1)-ones(LSIG,1)*mean(SIG(:,im1));
        signal2_centered=SIG(:,im2)-ones(LSIG,1)*mean(SIG(:,im2));
        Tfft_sec     = taustationary_sec/ratioDFT2SCP;
        Lfft         = fix(Tfft_sec*Fs_Hz);
        GREF         = ones(Lfft,1);
        [allSDs, time_sec, frqsFFT_Hz] = ...
            estimSCP(signal1_centered,signal2_centered,GREF,0.5, ...
            ratioDFT2SCP, 0, Fs_Hz, 'hann');
        allMSC{cp} = [allSDs.MSC];
    end
end
allMSCsort = allMSC(indsortdistance);

id1 = find(frqsFFT_Hz<=bandwidthMSC_Hz(1),1,'last');
id2 = find(frqsFFT_Hz<=bandwidthMSC_Hz(2),1,'last');

listindfreq  = (id1:id2);
frqsselected = frqsFFT_Hz(id1:id2);
Lf           = length(listindfreq);

%%
% addfig = 20;
% figure(1+addfig)
% 
% for ip=1:combi
%     subplot(7,7,ip)
%     pcolor(time_sec.SD/3600,frqsselected, allMSCsort{ip}(id1:id2,:))
%     shading flat
% end
%%
% slope = zeros(Lf,2);
% allMSC10 = zeros(combi,1);
% for ifreq=1:Lf
%     freq_ii = listindfreq(ifreq);
%     for ip=1:combi
%         allMSC10(ip) = (mean(allMSCsort{ip}(freq_ii,allMSCsort{ip}(freq_ii,:)>0.)));
%     end
%  
% %   plot(sortdistances,allMSC10,'x')
% %   title(sprintf('freq = %5.2f',frqsFFT_Hz(freq_ii)))
% %   pause
%     logallMSC = log(allMSC10);
%     HH        = [ones(length(sortdistances),1) sortdistances]; %
%     alphareg  = HH \ logallMSC;
%     slope(ifreq,:) = alphareg;%  ;
%     figure(2+addfig)
%     subplot(211)
%     plot(sortdistances, [HH*alphareg logallMSC],'.-')
%     hold on
%     subplot(212)
%     plot(sortdistances, abs(logallMSC-HH*alphareg) ./ abs(logallMSC))
%     hold on
%     pause
% end
%%
% clf
% LSCP = length(time_sec.SD);
% for ip=1:combi
%     aux = NaN(Lf,LSCP);
%     indaux          = find(allMSCsort{ip}(listindfreq,:)>0);
%     aux(indaux)     = allMSCsort{ip}(indaux);
%     logallMSCasfreq = log((aux));
%         plot(frqsFFT_Hz(listindfreq).^2, nanmean(logallMSCasfreq,2),'.-','color',allcolors(ip,1))
% hold on
% end
% % set(gca,'ylim',[log(0.3) log(0.8)])
% hold off
%%
clf
LSCP = length(time_sec.SD);
logaux = NaN(combi,LSCP,Lf);
for ifq=1:Lf
    aux = NaN(combi,LSCP);
    for ip=1:combi
        aux(ip,:) = allMSCsort{ip}(listindfreq(ifq),:);
    end
    logaux(:,:,ifq) = log(aux);
        plot(frqsFFT_Hz(listindfreq(ifq))^2*sortdistances .^2, ...
            nanmean(logaux(:,:,ifq),2),'.-','color',allcolors(ifq,1))
hold on
end
hold off

% set(gca,'ylim',[log(0.3) log(0.8)])
 set(gca,'xlim',[0 sortdistances(25)* sortdistances(25)*frqsFFT_Hz(listindfreq(ifq))^2])
 
%  imagesc(sortdistances .^2,frqsFFT_Hz(listindfreq) .^2,squeeze((nanmean(logaux,2))))
 
%%
% figure(2+addfig)
%     subplot(211)
% hold off
%     subplot(212)
% hold off
% figure(3+addfig)
% % subplot(211)
% plot(frqsFFT_Hz(listindfreq)',-slope(:,2)','.-');% .* frqsFFT_Hz(listindfreq)
% 
% [ones(Lf,1) frqsFFT_Hz(listindfreq)']\slope(:,2)
% grid on
% % subplot(212)
% % semilogx(frqsFFT_Hz(listindfreq)',slope(:,2)' ./ frqsFFT_Hz,'.-')
% % grid on
% % subplot(313)
% % semilogx(frqsFFT_Hz(listindfreq)',slope(:,3)','.-')
% % grid on
% %         figure(indexofSTA2)
% %         pcolor(time_sec.SD/60/60, frqsFFT_Hz(1:400),allMSC(1:400,:))
% %         shading flat
% %         set(gca,'yscale','log')
% %         colorbar
% 
