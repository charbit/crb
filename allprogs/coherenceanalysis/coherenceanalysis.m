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

addpath ../progsCRB/
%=================================================================

%=================================================================
stationnumber         = '31';

%=================================================================

%=================================================================
directorydatafromIDC  = sprintf('../../../../AAdataI%s/',stationnumber);

filenames              = dir(sprintf('%s*.mat',directorydatafromIDC));
nbfiles                = length(filenames);
Fs_Hz                  = 20;
Ts_sec                 = 1/Fs_Hz;
BWFilter4anal_Hz       = [0.01 4];
[forward,  reverse]    = butter(2,2*BWFilter4anal_Hz/Fs_Hz);
usefilterflag          = 1;
%========= MSC analysis
timeofanalysis_sec     = 500;
ratioDFT2SCP4average   = 5;
overlappingFFTrate     = 0.5;

for ifile = 1
    filename1_ii = filenames(ifile).name;filename1_ii
    cdload       = sprintf('load(''%s%s'');',directorydatafromIDC,filename1_ii);
    eval(cdload)
    
    xsensors_m   = observations.xsensors_m.coordinates;
    Msensors     = size(xsensors_m,1);
    xsensors_m   = xsensors_m-ones(Msensors,1)*xsensors_m(1,:);
    
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
    
    orientations = zeros(combi,1);
    cp           = 0;
   for i1=1:Msensors-1
        for i2=i1+1:Msensors
            cp=cp+1;
            orientations(cp) = atan(xsensors_m(i1,2)/xsensors_m(i2,2))*180/pi;
        end
    end
    sortorientations = orientations(indsortdistance);
    
    
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
    allSpectrum = cell(Msensors,1);
    cp=0;
    for im1 = 1:Msensors-1
        signal1_centered = SIGcentered(:,im1);
        for im2 = im1+1:Msensors
            cp = cp+1;
            signal2_centered = SIGcentered(:,im2);
            Tfft_sec     = timeofanalysis_sec/ratioDFT2SCP4average;
            Lfft         = fix(Tfft_sec*Fs_Hz);
            GREF         = ones(Lfft,1);
            [allSDs, time_sec, frqsFFT_Hz] = ...
                estimSCP(signal1_centered,signal2_centered,GREF, ...
                overlappingFFTrate, ...
                ratioDFT2SCP4average, 0, Fs_Hz, 'hann');
            allMSC{cp} = [allSDs.MSC];
        end
        allSpectrum{im1} = [allSDs.UU];
    end
    
    allSpectrum{Msensors} = [allSDs.RR];
    allMSCsort = allMSC(indsortdistance);
    LSCP = length(time_sec.SD);
    bandwidthdisplay_Hz   = [0.01 4];
    id1a         = find(frqsFFT_Hz<=bandwidthdisplay_Hz(1),1,'last');
    id2a         = find(frqsFFT_Hz<=bandwidthdisplay_Hz(2),1,'last');
    listindfreqa = (id1a:id2a);
    frqsselected_Hza   = frqsFFT_Hz(id1a:id2a);
    %%
    figure(101)
    im=1;
    subplot(2,1,1),
    plot(time_sec.signals/3600,SIGcentered(:,im),'color',0.5*[1 1 1]),
    set(gca,'xlim',[0 20]),
    set(gca,'ylim',[-10 10]),
    set(gca,'fontnam','times','fontsize',10)
    if or(im==9,im==10)
        xlabel('times - [H]')
    end
    title(sprintf('H%i : filtered in [%4.2f Hz - %4.2f Hz]',im, ...
        BWFilter4anal_Hz(1),BWFilter4anal_Hz(2)))
    %
    subplot(2,1,2),
    pcolor(time_sec.SD/3600,frqsselected_Hza,10*log10(allSpectrum{im}(listindfreqa,:))),
    shading flat
    set(gca,'ylim',[0 4]);
    set(gca,'xlim',[0 20]),
    set(gca,'yscale','log')
    set(gca,'fontnam','times','fontsize',10)
    xlabel('times - [H]')
    ylabel('freq. - [Hz]')
    title(sprintf('Spectral analysis in %i second time window',timeofanalysis_sec))
    
    
    HorizontalSize = 10;
    VerticalSize   = 8;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a4');
    %     set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    
    set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
    set(gcf, 'InvertHardCopy', 'off');
    
    figure(101)
    if BWFilter4anal_Hz(1)==0.1
        printname101 = sprintf('../../figures/tempspectanalysisH1HIGHI%s%s.eps',stationnumber,filename1_ii(1:8));
    else
        printname101 = sprintf('../../figures/tempspectanalysisH1LOWI%s%s.eps',stationnumber,filename1_ii(1:8));
    end
    fileprintepscmd101 = sprintf('print -depsc -loose %s',printname101);
    fileeps2pdf101 = sprintf('!epstopdf  %s',printname101);
    rmeps101 = sprintf('!rm %s',printname101);
%          eval(fileprintepscmd101)
%          eval(fileeps2pdf101)
%          eval(rmeps101)

%%
%
%     figure(100)
%     for im=1:Msensors,
%         subplot(Msensors/2,4,2*im-1),
%         plot(time_sec.signals/3600,SIGcentered(:,im),'color',0.5*[1 1 1]),
%         set(gca,'xlim',[0 20]),
%         set(gca,'ylim',[-3 3]),
%         set(gca,'fontnam','times','fontsize',10)
%         if or(im==9,im==10)
%             xlabel('times - [H]')
%         end
%         title(sprintf('H%i',im))
%         %
%         subplot(Msensors/2,4,2*im),
%         pcolor(time_sec.SD/3600,frqsselected_Hza,10*log10(allSpectrum{im}(listindfreqa,:))),
%         shading flat
%         set(gca,'ylim',[0 4]);%bandwidthdisplay_Hz),
%         set(gca,'xlim',[0 20]),
%         set(gca,'yscale','log')
%         set(gca,'fontnam','times','fontsize',10)
%         if or(im==9,im==10)
%             xlabel('times - [H]')
%         end
%         ylabel('freq. - [Hz]')
%     end
%     HorizontalSize = 28;
%     VerticalSize   = 26;
%     set(gcf,'units','centimeters');
%     set(gcf,'paperunits','centimeters');
%     set(gcf,'PaperType','a4');
%     %     set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
%     set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
%     
%     set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
%     set(gcf, 'InvertHardCopy', 'off');
%     
%     figure(100)
%     printname100 = sprintf('../../figures/tempspectanalysisI%s%s.eps',stationnumber,filename1_ii(1:8));
%     fileprintepscmd100 = sprintf('print -depsc -loose %s',printname100);
%     fileeps2pdf100 = sprintf('!epstopdf  %s',printname100);
%     rmeps100 = sprintf('!rm %s',printname100);
% %          eval(fileprintepscmd100)
% %          eval(fileeps2pdf100)
% %          eval(rmeps100)
    
    
    %%
    for HIGHBAND = [0,1]
        %===================================================================
        if HIGHBAND
            ifig = 22;
            bandwidthMSC_Hz   = [0.2 0.3];
            maxvarexplic_Hz2km2 = 0.15;
            nameprint = sprintf('../../figures/coherence2nearestI%s%sHIGH.eps',stationnumber,filename1_ii(1:8));
        else
            ifig = 23;
            bandwidthMSC_Hz   = [0.05 0.1];
            maxvarexplic_Hz2km2 = 0.02;
            nameprint = sprintf('../../figures/coherence2nearestI%s%sLOW.eps',stationnumber,filename1_ii(1:8));
        end
        
        idf1               = find(frqsFFT_Hz<=bandwidthMSC_Hz(1),1,'last');
        idf2               = find(frqsFFT_Hz<=bandwidthMSC_Hz(2),1,'last');
        listindfreq        = (idf1:idf2);
        frqsselected_Hz    = frqsFFT_Hz(idf1:idf2);
        nbfreq4MSC         = length(listindfreq);
        
        MSCtheresholdseed   = 0.9;
        %(bandwidthMSC_Hz(2)^2)*(sortdistances(combi)^2)*1e-6;
        %===================================================================
        figure(ifig)
        clf
        %=======
        subplot(421)
        plot((0:LSIG-1)/Fs_Hz/3600, SIG(:,sortdistances(1,2))/max(max(SIG)))
        set(gca, 'ylim',[-1.5 1.5])
        
        subplot(423)
        plot((0:LSIG-1)/Fs_Hz/3600, SIG(:,sortdistances(1,3))/max(max(SIG)))
        set(gca, 'ylim',[-1.5 1.5])
        
        
        
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
                figure(ifig)
                subplot(223)
                hold on
                plot(timeselect_samples{ifq}/Fs_Hz/3600,frqsFFT_Hz(listindfreq(ifq)),'.',...
                    'color',allcolors(ifqcolor,1))
                %=======
                
                explicativevar_Hz2km2     = frq_ifq_Hz2*(sortdistances(:,1) .^2)/1e6;
                explicativevarsave_Hz2km2 = [explicativevarsave_Hz2km2;explicativevar_Hz2km2];
                %=======
                
                subplot(122)
                semilogy(explicativevar_Hz2km2, exp(meanlogaux),'.-','color',allcolors(ifqcolor,1))
                hold on
                %                         semilogy(explicativevar_Hz2km2, exp(logaux(:,:,ifq)),'--','color',0.7*[1 1 1])
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
        figure(ifig)
        subplot(421)
        set(gca,'xlim',[0 20])
        set(gca,'ylim',0.5*[-1 1])
        set(gca,'ytick',[])
        grid on
        %     xlabel('time - [hour]','fontname','times','fontsize',10)
        set(gca,'fontname','times','fontsize',10)
        title(sprintf('day: %s/%s/%s',filename1_ii(1:4),filename1_ii(5:6),filename1_ii(7:8)))
        
        subplot(423)
        set(gca,'xlim',[0 20])
        set(gca,'ylim',0.5*[-1 1])
        set(gca,'ytick',[])
        grid on
        xlabel('time - [hour]','fontname','times','fontsize',10)
        set(gca,'fontname','times','fontsize',10)
        
        subplot(223)
        hold off
        set(gca,'xlim',[0 20])
        set(gca,'ylim',[bandwidthMSC_Hz(1)*0.9 bandwidthMSC_Hz(2)*1.1])
        set(gca,'yscale','lin')
        grid on
        ylabel('frequency - [Hz]','fontname','times','fontsize',10)
        set(gca,'fontname','times','fontsize',10)
        set(gca,'box','on')
        
        subplot(122)
        
        grid on
        %      set(gca,'xlim',[0 0.4]);%maxvarexplic_Hz2km2])
        set(gca,'ylim',[5e-2 1])
        set(gca,'fontname','times','fontsize',10)
        xlabel(sprintf('F%s x d%s - [Hz%s x km%s]',exp2,exp2,exp2,exp2),'fontname','times','fontsize',10)
        ylabel('MSC','fontname','times','fontsize',10)
        
        hold on
        semilogy(sortexplicativevarsave_Hz2km2([1 end]),...
            exp([ones(2,1) sortexplicativevarsave_Hz2km2([1 end])]*coeffs),'k','linew',1.5)
        hold off
        title(sprintf('LOC - decay = %4.2e Hz%s x m%s\nSNR = %4.2f dB',coeffs(2)*1e-6,exp2,exp2,-10*log10(exp(-coeffs(1)/2)-1)), ...
            'fontname','times','fontsize',10)
        %===== we print
        HorizontalSize = 20;
        VerticalSize   = 14;
        set(gcf,'units','centimeters');
        set(gcf,'paperunits','centimeters');
        set(gcf,'PaperType','a4');
        %     set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
        set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
        
        set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
        set(gcf, 'InvertHardCopy', 'off');
        
        figure(ifig)
        fileprintepscmd = sprintf('print -depsc -loose %s',nameprint);
        fileeps2pdf = sprintf('!epstopdf %s',nameprint);
        rmeps = sprintf('!rm  %s',nameprint);
        figure(ifig)
%         eval(fileprintepscmd)
%         eval(fileeps2pdf)
%         eval(rmeps)
    end
end