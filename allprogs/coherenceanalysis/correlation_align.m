%===== display the correlation wrt 1 for the full duration
% must be runned after the program coherenceanalysis.m
%===========================================================
Lalign_sec         = 500;
Lcorr_pts          = fix(Lalign_sec*Fs_Hz);
bandwidthMSC_Hz    = [0.2 0.3];
idf1               = find(frqsFFT_Hz<=bandwidthMSC_Hz(1),1,'last');
idf2               = find(frqsFFT_Hz<=bandwidthMSC_Hz(2),1,'last');
listindfreq        = (idf1:idf2);
frqsselected_Hz    = frqsFFT_Hz(idf1:idf2);
nbfreq4MSC         = length(listindfreq);
azimuth_deg        = NaN(LSCP,1);
horizvelocity_mps  = NaN(LSCP,1);
cormax             = zeros(Msensors-1,LSCP);
MSCsel             = zeros(Msensors-1,LSCP);
for iT=1:LSCP
    idt1 = fix(time_sec.SD(iT)*Fs_Hz)+1;
    idt2 = idt1+Lcorr_pts-1;
    
    signalsample = zeros(Lcorr_pts,Msensors);
    for im=1:Msensors
        signalsample(:,im) = SIGcentered(idt1:idt2,im);
    end
    [xalign, taupts, signal_notalign, cormax(:,iT)]   = ...
        alignmentwrt1(signalsample,1,Lcorr_pts,1,bandwidthMSC_Hz/Fs_Hz);

    incorrmax       = find(cormax(:,iT)>0.8);
    for im=1:Msensors-1
    MSCsel(im,iT)    = mean([allMSC{im}(idf1:idf2,iT)]);
    end
    if  length(incorrmax)>3
    thetaSN = xsensors_m(incorrmax+1,1:2)\(taupts(incorrmax+1)/Fs_Hz);
    horizvelocity_mps(iT) = 1/norm(thetaSN);
    azimuth_deg(iT) = atan2(thetaSN(2),thetaSN(1))*180/pi;
    end
end
%%
figure(1)
subplot(411); plot(time_sec.SD/3600, horizvelocity_mps)
subplot(412); plot(time_sec.SD/3600, azimuth_deg) 
% subplot(413); plot(time_sec.SD/3600, MSCsel) 
% subplot(414); plot(time_sec.SD/3600, cormax) 
subplot(212); plot(cormax, MSCsel,'.') 