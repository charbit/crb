%===== display the correlation wrt 1 for the full duration
% must be runned after the program coherenceanalysis.m
%===========================================================
Lalign_sec = 300;
bandwidthMSC_Hz   = [0.06 0.08];
Lcorr_pts = fix(Lalign_sec*Fs_Hz);
idf1               = find(frqsFFT_Hz<=bandwidthMSC_Hz(1),1,'last');
idf2               = find(frqsFFT_Hz<=bandwidthMSC_Hz(2),1,'last');
listindfreq        = (idf1:idf2);
frqsselected_Hz    = frqsFFT_Hz(idf1:idf2);
nbfreq4MSC         = length(listindfreq);
for iT=1:LSCP
    idt1 = fix(time_sec.SD(iT)*Fs_Hz)+1;
    idt2 = idt1+Lcorr_pts-1;
    
    signalsample = zeros(Lcorr_pts,Msensors);
    for im=1:Msensors
        signalsample(:,im) = SIGcentered(idt1:idt2,im);
    end
    
    
    
    [xalign, taupts, signal_notalign, cormax]   = ...
        alignmentwrt1(signalsample,1,Lcorr_pts,bandwidthMSC_Hz/Fs_Hz);
    
    %===============
    figure(1000)
    subplot(311);
    plot((1:Lalign_sec*Fs_Hz)/Fs_Hz/60,signal_notalign)
    %     set(gca,'xlim',[6000 7000])
    subplot(312);
    plot((1:Lalign_sec*Fs_Hz)/Fs_Hz/60,xalign(:,:))
    %     set(gca,'xlim',[6000 7000])
    drawnow
    incorrmax = find(cormax>0.8);
    if length(incorrmax)>2
        thetaSN = xsensors_m(incorrmax+1,1:2)\(taupts(incorrmax+1)/Fs_Hz);
        1/norm(thetaSN)
        atan2(thetaSN(2),thetaSN(1))*180/pi
        subplot(313);
        plot(distances(1:Msensors-1,1),cormax,'o')
        hold on
    else
        disp('NAN')
    end
end
subplot(313);
hold off