%====================== coherenceanalysis.m =============================
% file contains OBSERVATIONS structure
% with 
%           .xsensors_m.coordinates : sensor coordinates in m
%           .xsensors_m.name
%           .data : data
% used functions: 
%   estimSCP.m
%======================
clear
exp2  = '^{-2}';
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.'; ...
    'bo';'ro';'mo';'co';'go';'ko';
    'bx';'rx';'mx';'cx';'gx';'kx'];

switch computer
    case 'GLNXA64'
        addpath /dvlscratch/SHI/users/charbit/ProjectIMS2015b/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
    otherwise
        addpath /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
end

%=================================================================

%=================================================================
stationnumber         = 31;

%=================================================================

%=================================================================
directorydatafromIDC  = sprintf('../../../../AAdataI%i/',stationnumber);

filenames              = dir(sprintf('%s*.mat',directorydatafromIDC));
nbfiles                = length(filenames);
Fs_Hz                  = 20;
Ts_sec                 = 1/Fs_Hz;
BWFilter4anal_Hz       = [0.01 5];
[forward,  reverse]    = butter(2,2*BWFilter4anal_Hz/Fs_Hz);
usefilterflag          = 1;
%========= MSC analysis
timeofanalysis_sec     = 500;
ratioDFT2SCP4average   = 5;
overlappingFFTrate     = 0.5;

for ifile = 10
    filename1_ii = filenames(ifile).name;filename1_ii
    cdload       = sprintf('load(''%s%s'');',directorydatafromIDC,filename1_ii);
    eval(cdload)
    
    xsensors_m   = observations.xsensors_m.coordinates;
    Msensors     = size(xsensors_m,1);
    cp           = 0;
    combi        = Msensors*(Msensors-1)/2;
    distances    = zeros(combi,3);
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
    %         orientations(cp) = atan((xsensors_m(i1,2)-xsensors_m(i2,2)) / 
    %         (xsensors_m(i1,2)-xsensors_m(i2,2));
    %     end
    % end
    % [sortdistances, indsortdistance] = sort(distances);
    
    
    LSIG = size(observations.data,1);
    SIG  = zeros(LSIG,Msensors);
    for im=1:Msensors
        SIG(:,im) = observations.data(:,im);
        if usefilterflag
            SIG(:,im) = filter(forward,reverse,SIG(:,im));
        end
    end
 
    SIGcentered = SIG - ones(LSIG,1)*mean(SIG,1);
    
    %=====================
    %======== MSC analysis
    %=====================
    allMSC = cell(combi,1);

    cp=0;
    for im1 = 1:Msensors-1
        for im2 = im1+1:Msensors
            cp = cp+1;
            signal1_centered=SIGcentered(:,im1);
            signal2_centered=SIGcentered(:,im2);
            Tfft_sec     = timeofanalysis_sec/ratioDFT2SCP4average;
            Lfft         = fix(Tfft_sec*Fs_Hz);
            GREF         = ones(Lfft,1);
            [allSDs, time_sec, frqsFFT_Hz] = ...
                estimSCP(signal2_centered,signal1_centered,GREF, ...
                overlappingFFTrate, ...
                ratioDFT2SCP4average, 0, Fs_Hz, 'hann');
            allMSC{cp} = [allSDs.MSC];
        end
    end
    allMSCsort = allMSC(indsortdistance);
        
    %%
    bandwidthMSC_Hz   = [0.08 0.12];
    id1               = find(frqsFFT_Hz<=bandwidthMSC_Hz(1),1,'last');
    id2               = find(frqsFFT_Hz<=bandwidthMSC_Hz(2),1,'last');
    listindfreq       = (id1:id2);
    frqsselected_Hz   = frqsFFT_Hz(id1:id2);
    nbfreq4MSC        = length(listindfreq);
    
    MSCtheresholdseed   = 0.9;
    maxvarexplic_Hz2km2 = 0.02;%(bandwidthMSC_Hz(2)^2)*(sortdistances(combi)^2)*1e-6;
    %===================================================================
    figure(ifile)
    clf
    %=======
    subplot(221)
    plot((0:LSIG-1)/Fs_Hz/3600, SIG(:,sortdistances(1,2:3))/max(max(SIG)))
    set(gca, 'ylim',[-1.5 1.5])
    
    LSCP = length(time_sec.SD);
    logaux = NaN(combi,LSCP,nbfreq4MSC);
    timeselect_samples = cell(nbfreq4MSC,1);
    meanlogauxsave = [];
    explicativevarsave_Hz2km2 = [];
    for ifq=1:nbfreq4MSC
        ifqcolor   = mod(ifq,length(allcolors))+1;
        frq_ifq_Hz = frqsselected_Hz(ifq);
        frq_ifq_Hz2 = frq_ifq_Hz^2;
        aux        = NaN(combi,LSCP);
        indselect = find(and(and(...
            allMSCsort{1}(listindfreq(ifq),:)>MSCtheresholdseed,...
            allMSCsort{2}(listindfreq(ifq),:)>0),...
            allMSCsort{3}(listindfreq(ifq),:)>0));
        if length(indselect)>1
            timeselect_samples{ifq} = time_sec.SD(indselect)*Fs_Hz;
            for ip=1:combi
                aux(ip,indselect) = allMSCsort{ip}(listindfreq(ifq),indselect);
            end
            logaux(:,:,ifq) = log(aux);
            meanlogaux =  nanmean(logaux(:,:,ifq),2);
            stdlogaux =  nanstd(logaux(:,:,ifq),[],2);
            meanlogauxsave = [meanlogauxsave;meanlogaux];
            
            %=======
            figure(ifile)
            subplot(223)
            hold on           
            plot(timeselect_samples{ifq}/Fs_Hz/3600,frqsFFT_Hz(listindfreq(ifq)),'.',...
                'color',allcolors(ifqcolor,1))

            explicativevar_Hz2km2     = frq_ifq_Hz2*(sortdistances(:,1) .^2)/1e6;
            explicativevarsave_Hz2km2 = [explicativevarsave_Hz2km2;explicativevar_Hz2km2];
            %=======
            figure(ifile)
            subplot(122)
            semilogy(explicativevar_Hz2km2, exp(meanlogaux),'.-','color',allcolors(ifqcolor,1))
            hold on
%             semilogy(explicativevar_Hz2km2, exp(logaux(:,:,ifq)),'--','color',0.7*[1 1 1])
%             semilogy(explicativevar_Hz2km2, exp(meanlogaux+stdlogaux),'--','color',0.3*[1 1 1])
%             semilogy(explicativevar_Hz2km2, exp(meanlogaux-stdlogaux),'--','color',0.3*[1 1 1])
        end
    end
    
    %=== regression on coefs
    [sortexplicativevarsave_Hz2km2, indsorteva] = ...
        sort(explicativevarsave_Hz2km2);
    sortmeanlogauxsave                          = meanlogauxsave(indsorteva);
    auxnum                                      = ...
        find(sortexplicativevarsave_Hz2km2<maxvarexplic_Hz2km2,1,'last');
    HH          = [ones(auxnum,1) sortexplicativevarsave_Hz2km2(1:auxnum)];
    coeffs      = HH\sortmeanlogauxsave(1:auxnum);
     
    %======================================================
    figure(ifile)
    subplot(221)
    set(gca,'xlim',[0 24])
    grid on
    xlabel('time - [hour]','fontname','times','fontsize',10)
        set(gca,'fontname','times','fontsize',10)

    subplot(223)
    hold off
    set(gca,'xlim',[0 24])
    set(gca,'ylim',[0.01 bandwidthMSC_Hz(2)*1.1])
    set(gca,'yscale','log')
    grid on
    ylabel('frequency - [Hz]','fontname','times','fontsize',10)
    set(gca,'fontname','times','fontsize',10)
    set(gca,'box','on')
    
    subplot(122)
    grid on
%     set(gca,'xlim',[0 maxvarexplic_Hz2km2])
    set(gca,'ylim',[1e-1 1])    
    set(gca,'fontname','times','fontsize',10)
    xlabel(sprintf('F2 x d2 - [Hz%s x km%s]',exp2,exp2),'fontname','times','fontsize',10)
    ylabel('MSC','fontname','times','fontsize',10)

    figure(ifile)
    subplot(122)
    hold on
    semilogy(sortexplicativevarsave_Hz2km2([1 end]),...
        exp([ones(2,1) sortexplicativevarsave_Hz2km2([1 end])]*coeffs),'k','linew',1.5)
    hold off
    title(sprintf('LOC - decay = %4.2e Hz%s x m%s',coeffs(2)*1e-6,exp2,exp2), ...
        'fontname','times','fontsize',10)    

    HorizontalSize = 24;
    VerticalSize   = 10;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a4');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    
    set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
    set(gcf, 'InvertHardCopy', 'off');
    fileprintepscmd = sprintf('print -depsc -loose ../../figures/coherenceI%i%s.eps',stationnumber,filename1_ii(1:8));
%     eval(fileprintepscmd)

end