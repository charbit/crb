%====================== estimationwithFB.m =============================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

switch computer
    case 'GLNXA64'
        addpath /dvlscratch/SHI/users/charbit/ProjectIMS2015b/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
    otherwise
        addpath /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
end

stationnumber         = 37;
directorydatafromIDC  = sprintf('../../../../AAdataI%i/',stationnumber);

filenames              = dir(sprintf('%s*.mat',directorydatafromIDC));
nbfiles                = length(filenames);
Fs_Hz                  = 20;
Ts_sec                 = 1/Fs_Hz;

for ifile = 16%:nbfiles
    filename1_ii = filenames(ifile).name;
    cdload = sprintf('load(''%s%s'');',directorydatafromIDC,filename1_ii);
    eval(cdload)
    
    xsensors_m = observations.xsensors_m.coordinates;
    Msensors   = size(xsensors_m,1);
    cp = 0;
    combi = Msensors*(Msensors-1)/2;
    distances = zeros(combi,3);
    for i1=1:Msensors-1
        for i2=i1+1:Msensors
            cp=cp+1;
            distances(cp,:) = [...
                norm(xsensors_m(i1,:)-xsensors_m(i2,:)) ...
                i1,i2];
        end
    end
    [sortdistances, indsortdistance] = sort(distances(:,1));
    sortdistances(:,2:3) = distances(indsortdistance,2:3);
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
    bandwidthFilter_Hz  = [0.02 1];
    [forward,  reverse] = butter(2,2*bandwidthFilter_Hz/Fs_Hz);
    for im = 1:Msensors
        SIG(:,im)       = filter(forward,reverse,SIG(:,im));
    end
    SIGcentered = SIG - ones(LSIG,1)*mean(SIG,1);
    % s1=SIGcentered(115000+(-10000:10000),1);
    % s2=SIGcentered(115000+(-10000:10000),2);
    % subplot(212)
    % mscohere(s1,s2,[],[],[],20)
    % set(gca,'xlim',[0 0.4])
    % subplot(211)
    % plot((0:20000)/Fs_Hz,[s1 s2])
    %
    %
    % return
    %
    %
    %     bandwidthMSC_Hz   = [0.05 0.2];
    taustationary_sec = 500 ;
    ratioDFT2SCP      = 5;
    allMSC = cell(combi,1);
    
    %=====================
    %=====================
    %=====================
    
    
    
    cp=0;
    for im1 = 1:Msensors-1
        for im2 = im1+1:Msensors
            cp = cp+1;
            signal1_centered=SIGcentered(:,im1);
            signal2_centered=SIGcentered(:,im2);
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
    
    cdsave = sprintf('save %s%s allMSCsort taustationary_sec ratioDFT2SCP SIGcentered sortdistances indsortdistance bandwidthFilter_Hz time_sec frqsFFT_Hz xsensors_m',directorydatafromIDC,filename1_ii);
   cdsave = sprintf('save %s%s',directorydatafromIDC,filename1_ii);
    eval(cdsave)
    
    %%
    bandwidthMSC_Hz   = [0.08 0.12];
    
    id1 = find(frqsFFT_Hz<=bandwidthMSC_Hz(1),1,'last');
    id2 = find(frqsFFT_Hz<=bandwidthMSC_Hz(2),1,'last');
    
    listindfreq  = (id1:id2);
    frqsselected = frqsFFT_Hz(id1:id2);
    Lf           = length(listindfreq);
    
    
    
    figure(ifile)
    clf
    subplot(311)
    plot((0:LSIG-1)/Fs_Hz/3600, SIG(:,sortdistances(1,2:3))/max(max(SIG)))
    set(gca, 'ylim',[-1.5 1.5])
    
    LSCP = length(time_sec.SD);
    logaux = NaN(combi,LSCP,Lf);
    timeselect_samples = cell(Lf,1);
    llll=linspace(0.8,1.2,Lf);
    for ifq=1:Lf
        aux = NaN(combi,LSCP);
        indselect = find(and(and(...
            allMSCsort{1}(listindfreq(ifq),:)>0.8,...
            allMSCsort{2}(listindfreq(ifq),:)>0.8),...
            allMSCsort{3}(listindfreq(ifq),:)>0.8));
        if length(indselect)>30
            timeselect_samples{ifq} = time_sec.SD(indselect)*Fs_Hz;
            for ip=1:combi
                aux(ip,indselect) = allMSCsort{ip}(listindfreq(ifq),indselect);
            end
            logaux(:,:,ifq) = log(aux);
            meanlogaux =  nanmean(logaux(:,:,ifq),2);
            stdlogaux =  nanstd(logaux(:,:,ifq),[],2);
            
            figure(ifile)
            subplot(312)
            hold on
            %         plot(timeselect_samples{ifq}/Fs_Hz/3600,llll(ifq)*ones(length(indselect),1),'.')
            
            plot(timeselect_samples{ifq}/Fs_Hz/3600,frqsFFT_Hz(listindfreq(ifq)),'.',...
                'color',allcolors(ifq,1))
            hold off
            
            figure(ifile)
            subplot(313)
            
            plot(frqsFFT_Hz(listindfreq(ifq))^2*sortdistances(:,1) .^2, ...
                meanlogaux,'.-','color',allcolors(ifq,1))
            hold on
                    plot(frqsFFT_Hz(listindfreq(ifq))^2*sortdistances(:,1) .^2, ...
                       meanlogaux+stdlogaux,'--','color',allcolors(ifq,1))
                    plot(frqsFFT_Hz(listindfreq(ifq))^2*sortdistances(:,1) .^2, ...
                       meanlogaux-+stdlogaux,'--','color',allcolors(ifq,1))
        end
    end
    figure(ifile)
    subplot(311)
    subplot(312)
    set(gca,'ylim',[0 0.5])
    subplot(313)
    hold off
    grid on
    set(gca,'ylim',[-2 0])
    set(gca,'xlim',[0 2e4])
    %    set(gca,'xlim',[0 sortdistances(25,1)* sortdistances(25,1)*frqsFFT_Hz(listindfreq(ifq))^2])
    
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
    
    %  imagesc(sortdistances .^2,frqsFFT_Hz(listindfreq) .^2,squeeze((nanmean(logaux,2))))
    
    %%1
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
end




% %%
% for ifile=1:nbfiles
%     figure(ifile)
%     xlabel('distances - [m]')
%     ylabel('log MSC')
%     hfig=20;
%     vfig=13;
%     set(gcf,'units','centimeters');
%     set(gcf,'paperunits','centimeters');
%     set(gcf,'PaperType','a4');
% %     set(gcf,'position',[8 10 hfig vfig]);
%     set(gcf,'paperposition',[0 0 hfig vfig]);
% end