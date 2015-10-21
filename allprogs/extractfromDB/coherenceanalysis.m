%====================== estimationwithFB.m =============================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];


% addpath /dvlscratch/SHI/users/charbit/ProjectIMS2015b/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
addpath /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
directorydatafromIDC  = '../AAdataI37/';

%==========================================================================
% the data are in a file with the name built as:
%                  sta1_Y2015_D239.mat
% meaning station number 1, year 2015, days 239 and 240
% these data have been extracted from IDC testbed
%==========================================================================
filenames              = dir(sprintf('%ssta*.mat',directorydatafromIDC));
nbfiles                = length(filenames);
Fs_Hz                  = 20;%records{ir}.Fs_Hz;
Ts_sec                 = 1/Fs_Hz;

%=====================
for indexofSTA1 = 1%:nbfiles-1
    
    for indexofSTA2 = 3%:nbfiles%indexofSTA1+1:2%nbfiles
        
        filename1_ii = filenames(indexofSTA1).name;
        cdload1 = sprintf('sig1 = load(''%s%s'');',directorydatafromIDC,filename1_ii);
        eval(cdload1)
        Lrecords1   = length(sig1.records);
        LL1_max=0;
        for ir =1:Lrecords1
            if length(sig1.records{ir}.data)>LL1_max
                ir1save = ir;
                LL1_max=length(sig1.records{ir}.data);
            end
        end
        signal1 = [sig1.records{ir1save}.data];
        
        
        filename2_ii = filenames(indexofSTA2).name;
        cdload2 = sprintf('sig2 = load(''%s%s'');',directorydatafromIDC,filename2_ii);
        eval(cdload2)
        Lrecords2   = length(sig2.records);
        LL2_max=0;
        for ir =1:Lrecords2
            if length(sig2.records{ir}.data)>LL2_max
                ir2save = ir;
                LL2_max=length(sig2.records{ir}.data);
            end
        end
        signal2 = [sig2.records{ir2save}.data];
        
        st1 = sig1.records{ir1save}.stime;
        st2 = sig2.records{ir2save}.stime;
        et1 = sig1.records{ir1save}.etime;
        et2 = sig2.records{ir2save}.etime;
        
        if and(st1>=st2,et1>=et2)
            ids2 = fix((st1-st2)*Fs_Hz)+1;
            ide1 = fix((et1-et2)*Fs_Hz)+1;
            signal1 = signal1(1:end-ide1+1);
            signal2 = signal2(ids2:end);
        end
        if and(st1>=st2,et1<et2)
            ids2 = fix((st1-st2)*Fs_Hz)+1;
            ide2 = fix((et2-et1)*Fs_Hz)+1;
            signal2 = signal2(ids2:end-ide2+1);
        end
        if and(st1<st2,et1>=et2)
            ids1 = fix((st2-st1)*Fs_Hz)+1;
            ide1 = fix((et1-et2)*Fs_Hz)+1;
            signal1 = signal1(ids1:end-ide1+1);
        end
        if and(st1<st2,et1<et2)
            ids1 = fix((st2-st1)*Fs_Hz)+1;
            ide2 = fix((et2-et1)*Fs_Hz)+1;
            signal1 = signal2(ids1:end);
            signal2 = signal2(1:end-ide2+1);
        end
        signal1_centered=signal1-ones(size(signal1,1),1)*mean(signal1);
        signal2_centered=signal2-ones(size(signal2,1),1)*mean(signal2);
        
        taustationary_sec = 400 ;
        ratioDFT2SCP = 5;
        Tfft_sec  = taustationary_sec/ratioDFT2SCP;
        Lfft      = fix(Tfft_sec*Fs_Hz);
        
        GREF = ones(Lfft,1);
        
        [allSDs, time_sec, frqsFFT_Hz] = ...
            estimSCP(signal1_centered,signal2_centered,GREF,0.5, ...
            ratioDFT2SCP, 0, Fs_Hz, 'hann');
        allMSC = [allSDs.MSC];
        figure(indexofSTA2)
        pcolor(time_sec.SD/60/60, frqsFFT_Hz(1:400),allMSC(1:400,:))
        shading flat
        set(gca,'yscale','log')
        colorbar
    end
end

