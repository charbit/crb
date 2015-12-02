% allMSCsort
% fixed frequency
selectedfrequency_Hz = 0.27;
indselectedfrequency = find(frqsFFT_Hz<=selectedfrequency_Hz,1,'last');

        indselect_t = find(and(and(...
            allMSCsort{1}(indselectedfrequency,:)>0.8,...
            allMSCsort{2}(indselectedfrequency,:)>0.7),...
            allMSCsort{3}(indselectedfrequency,:)>0.6));

figure(1)

for ip=4:combi
    aux = NaN(LSCP,1);
    aux(indselect_t) = allMSCsort{ip}(indselectedfrequency,indselect_t)';
    plot(time_sec.SD/3600,aux,'.-','color',allcolors(ip,1))
    hold on
end
hold off


 